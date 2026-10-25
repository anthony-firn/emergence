//! GPU port of the IMEX ARS(3,4,3) step with implicit Brinkman penalization
//! (math and notation in `spectral_core::imex`).
//!
//! Memory: instead of storing every stage's N and L, each future stage's
//! right-hand side and the final combination are accumulated as soon as a
//! stage is finished, in accumulators held by the IMEX workspace.
//! That needs only one extra set each for N_i, L_i and a spectral scratch.

use super::*;

pub(crate) const CG_TOL: f64 = 2e-6;
pub(crate) const CG_MAX_ITERS: usize = 400;

pub(crate) struct ImexKernels {
    accum4: CudaFunction,
    hinv: CudaFunction,
    hinv_inplace: CudaFunction,
    sub_hinv: CudaFunction,
    precond: CudaFunction,
    stage_lin: CudaFunction,
    mask_sub: CudaFunction,
    apply_finish: CudaFunction,
    mask_zero: CudaFunction,
    axpy: CudaFunction,
    xpby: CudaFunction,
    sub: CudaFunction,
    dot3: CudaFunction,
    sum3: CudaFunction,
}

impl ImexKernels {
    pub(crate) fn load(m: &Arc<cudarc::driver::CudaModule>) -> Result<Self> {
        Ok(Self {
            accum4: m.load_function("accum4")?,
            hinv: m.load_function("hinv")?,
            hinv_inplace: m.load_function("hinv_inplace")?,
            sub_hinv: m.load_function("sub_hinv")?,
            precond: m.load_function("precond_combine")?,
            stage_lin: m.load_function("stage_lin")?,
            mask_sub: m.load_function("mask_sub")?,
            apply_finish: m.load_function("apply_finish")?,
            mask_zero: m.load_function("mask_zero")?,
            axpy: m.load_function("axpy_f")?,
            xpby: m.load_function("xpby_f")?,
            sub: m.load_function("sub_f")?,
            dot3: m.load_function("dot3_partial")?,
            sum3: m.load_function("sum3_partial")?,
        })
    }
}

pub(crate) struct ImexDev {
    nl: [CudaSlice<float2>; 3],
    lin: [CudaSlice<float2>; 3],
    spec2: [CudaSlice<float2>; 3],
    b: [CudaSlice<f32>; 3],
    r: [CudaSlice<f32>; 3],
    z: [CudaSlice<f32>; 3],
    p: [CudaSlice<f32>; 3],
    ap: [CudaSlice<f32>; 3],
    g_prev: [[CudaSlice<f32>; 3]; 3],
    /// accumulators for the right-hand sides of stages 2..4 and uⁿ⁺¹
    racc: [[CudaSlice<float2>; 3]; 3],
    fin: [CudaSlice<float2>; 3],
    pub(crate) last_iters: [usize; 3],
    pub(crate) last_resid: [f64; 3],
    pub(crate) last_force: [f64; 3],
}

impl ImexDev {
    fn new(s: &Arc<CudaStream>, ns: usize, nr: usize) -> Result<Self> {
        Ok(Self {
            nl: alloc3(s, ns)?,
            lin: alloc3(s, ns)?,
            spec2: alloc3(s, ns)?,
            b: alloc3(s, nr)?,
            r: alloc3(s, nr)?,
            z: alloc3(s, nr)?,
            p: alloc3(s, nr)?,
            ap: alloc3(s, nr)?,
            g_prev: [alloc3(s, nr)?, alloc3(s, nr)?, alloc3(s, nr)?],
            racc: [alloc3(s, ns)?, alloc3(s, ns)?, alloc3(s, ns)?],
            fin: alloc3(s, ns)?,
            last_iters: [0; 3],
            last_resid: [0.0; 3],
            last_force: [0.0; 3],
        })
    }
}

fn accum(
    w: &Work,
    ik: &ImexKernels,
    out: &mut [CudaSlice<float2>; 3],
    xs: [&[CudaSlice<float2>; 3]; 2],
    c: [f64; 2],
) -> Result<()> {
    let ns = w.ns() as i64;
    let (c0, c1) = (c[0] as f32, c[1] as f32);
    let z = 0.0f32;
    for comp in 0..3 {
        unsafe {
            w.stream
                .launch_builder(&ik.accum4)
                .arg(&mut out[comp])
                .arg(&xs[0][comp]).arg(&xs[1][comp]).arg(&xs[0][comp]).arg(&xs[0][comp])
                .arg(&c0).arg(&c1).arg(&z).arg(&z)
                .arg(&ns)
                .launch(cfg(ns as usize))?;
        }
    }
    Ok(())
}

fn dot3(w: &mut Work, ik: &ImexKernels, a: &[CudaSlice<f32>; 3], b: &[CudaSlice<f32>; 3]) -> Result<f64> {
    let nr = w.nr() as i64;
    let lc = LaunchConfig {
        grid_dim: (REDUCE_BLOCKS, 1, 1),
        block_dim: (BLOCK, 1, 1),
        shared_mem_bytes: BLOCK * 8,
    };
    unsafe {
        w.stream
            .launch_builder(&ik.dot3)
            .arg(&a[0]).arg(&a[1]).arg(&a[2])
            .arg(&b[0]).arg(&b[1]).arg(&b[2])
            .arg(&mut w.reduce_f64)
            .arg(&nr)
            .launch(lc)?;
    }
    let parts = w.stream.clone_dtoh(&w.reduce_f64)?;
    Ok(parts[..REDUCE_BLOCKS as usize].iter().sum())
}

/// Solve (I − aL)U = r on the device. Returns (iterations, residual, force).
fn implicit_solve(
    w: &mut Work,
    ik: &ImexKernels,
    dev: &mut ImexDev,
    nu: f64,
    r: &[CudaSlice<float2>; 3],
    a: f64,
    out: &mut [CudaSlice<float2>; 3],
    stage: usize,
    keep_mean: bool,
    surf: Option<(&mut super::surface_gpu::DevSurface, &super::surface_gpu::NufftKernels)>,
) -> Result<(usize, f64, [f64; 3])> {
    let (n, nh) = (w.n, w.nh);
    let (ns, nr) = (w.ns(), w.nr());
    let nr64 = nr as i64;
    let anu = (a * nu) as f32;
    let s = w.stream.clone();
    for c in 0..3 {
        unsafe {
            s.launch_builder(&ik.hinv).arg(&mut out[c]).arg(&r[c]).arg(&anu).arg(&n).arg(&nh).launch(cfg(ns))?;
        }
    }
    if let Some((sd, nk)) = surf {
        return if sd.ext.is_some() {
            super::surface_gpu::extension_solve(w, nk, sd, nu, a, out, &mut dev.nl, &mut dev.lin, stage, keep_mean)
        } else {
            super::surface_gpu::surface_solve(w, nk, sd, nu, a, out, &mut dev.nl, &mut dev.lin, stage, keep_mean)
        };
    }
    if !w.penal.as_ref().is_some_and(|p| p.implicit) {
        return Ok((0, 0.0, [0.0; 3]));
    }
    // Take the penalization out so `w` can be borrowed mutably below; it is
    // put back before returning (an error here is fatal for the solver).
    let pen = w.penal.take().unwrap();
    let (eta, us) = (pen.eta, pen.u_solid);
    let chi = &pen.chi;
    let scale = 1.0f32 / nr as f32;
    let (eta_a, dbar) = ((eta / a) as f32, (eta / a) as f32);
    let zero_i = 0i32;
    let keep: i32 = keep_mean as i32;
    let zf = 0.0f32;

    // b = χ>0 · (F⁻¹H⁻¹r̂ − u_s)
    for c in 0..3 {
        Work::to_phys(&s, &w.c2r, &out[c], &mut w.tmp[c], &mut dev.b[c])?;
        let usc = us[c];
        unsafe {
            s.launch_builder(&ik.mask_sub).arg(&mut dev.b[c]).arg(chi).arg(&usc).arg(&nr64).launch(cfg(nr))?;
        }
    }

    // A x = D x + χ>0 · F⁻¹ H⁻¹ P F x, written into `ax`
    macro_rules! apply {
        ($x:expr, $ax:expr) => {{
            for c in 0..3 {
                w.r2c.exec_r2c(&$x[c], &mut w.tmp[c])?;
            }
            {
                let [t0, t1, t2] = &mut w.tmp;
                unsafe {
                    s.launch_builder(&w.k.project)
                        .arg(&mut *t0).arg(&mut *t1).arg(&mut *t2)
                        .arg(&dev.spec2[0]).arg(&dev.spec2[1]).arg(&dev.spec2[2])
                        .arg(&zero_i).arg(&scale).arg(&n).arg(&nh)
                        .arg(&keep).arg(&zf).arg(&zf).arg(&zf)
                        .arg(&w.kc2).arg(&zf).arg(&zf).arg(&zf)
                        .launch(cfg(ns))?;
                }
            }
            for c in 0..3 {
                unsafe {
                    s.launch_builder(&ik.hinv_inplace).arg(&mut w.tmp[c]).arg(&anu).arg(&n).arg(&nh).launch(cfg(ns))?;
                }
                w.c2r.exec_c2r(&mut w.tmp[c], &mut $ax[c])?;
                unsafe {
                    s.launch_builder(&ik.apply_finish)
                        .arg(&mut $ax[c]).arg(&$x[c]).arg(chi).arg(&eta_a).arg(&nr64)
                        .launch(cfg(nr))?;
                }
            }
        }};
    }
    // z = χ>0 · F⁻¹ M F r   (restricted Fourier inverse of D̄ + S)
    macro_rules! precond {
        ($rr:expr, $zz:expr) => {{
            for c in 0..3 {
                w.r2c.exec_r2c(&$rr[c], &mut w.tmp[c])?;
                s.memcpy_dtod(&w.tmp[c], &mut dev.spec2[c])?;
            }
            {
                let [t0, t1, t2] = &mut dev.spec2;
                unsafe {
                    s.launch_builder(&w.k.project)
                        .arg(&mut *t0).arg(&mut *t1).arg(&mut *t2)
                        .arg(&w.tmp[0]).arg(&w.tmp[1]).arg(&w.tmp[2])
                        .arg(&zero_i).arg(&scale).arg(&n).arg(&nh)
                        .arg(&keep).arg(&zf).arg(&zf).arg(&zf)
                        .arg(&w.kc2).arg(&zf).arg(&zf).arg(&zf)
                        .launch(cfg(ns))?;
                }
            }
            for c in 0..3 {
                unsafe {
                    s.launch_builder(&ik.precond)
                        .arg(&mut w.tmp[c]).arg(&dev.spec2[c]).arg(&anu).arg(&dbar).arg(&scale).arg(&n).arg(&nh)
                        .launch(cfg(ns))?;
                }
                w.c2r.exec_c2r(&mut w.tmp[c], &mut $zz[c])?;
                unsafe {
                    s.launch_builder(&ik.mask_zero).arg(&mut $zz[c]).arg(chi).arg(&nr64).launch(cfg(nr))?;
                }
            }
        }};
    }

    let mut g = std::mem::replace(&mut dev.g_prev[stage], alloc3(&s, 1)?);
    apply!(g, dev.ap);
    for c in 0..3 {
        unsafe {
            s.launch_builder(&ik.sub).arg(&mut dev.r[c]).arg(&dev.b[c]).arg(&dev.ap[c]).arg(&nr64).launch(cfg(nr))?;
        }
    }
    let bnorm = dot3(w, ik, &dev.b, &dev.b)?.sqrt().max(1e-30);
    precond!(dev.r, dev.z);
    for c in 0..3 {
        s.memcpy_dtod(&dev.z[c], &mut dev.p[c])?;
    }
    let mut rz = dot3(w, ik, &dev.r, &dev.z)?;
    let mut iters = 0;
    let mut rnorm = dot3(w, ik, &dev.r, &dev.r)?.sqrt();
    while rnorm > CG_TOL * bnorm && iters < CG_MAX_ITERS {
        apply!(dev.p, dev.ap);
        let alpha = (rz / dot3(w, ik, &dev.p, &dev.ap)?) as f32;
        let malpha = -alpha;
        for c in 0..3 {
            unsafe {
                s.launch_builder(&ik.axpy).arg(&mut g[c]).arg(&dev.p[c]).arg(&alpha).arg(&nr64).launch(cfg(nr))?;
                s.launch_builder(&ik.axpy).arg(&mut dev.r[c]).arg(&dev.ap[c]).arg(&malpha).arg(&nr64).launch(cfg(nr))?;
            }
        }
        precond!(dev.r, dev.z);
        let rz_new = dot3(w, ik, &dev.r, &dev.z)?;
        let beta = (rz_new / rz) as f32;
        rz = rz_new;
        for c in 0..3 {
            unsafe {
                s.launch_builder(&ik.xpby).arg(&mut dev.p[c]).arg(&dev.z[c]).arg(&beta).arg(&nr64).launch(cfg(nr))?;
            }
        }
        rnorm = dot3(w, ik, &dev.r, &dev.r)?.sqrt();
        iters += 1;
    }

    // Û = H⁻¹ r̂ − H⁻¹ P ĝ
    for c in 0..3 {
        w.r2c.exec_r2c(&g[c], &mut w.tmp[c])?;
    }
    {
        let [t0, t1, t2] = &mut w.tmp;
        unsafe {
            s.launch_builder(&w.k.project)
                .arg(&mut *t0).arg(&mut *t1).arg(&mut *t2)
                .arg(&dev.spec2[0]).arg(&dev.spec2[1]).arg(&dev.spec2[2])
                .arg(&zero_i).arg(&scale).arg(&n).arg(&nh)
                        .arg(&keep).arg(&zf).arg(&zf).arg(&zf)
                        .arg(&w.kc2).arg(&zf).arg(&zf).arg(&zf)
                .launch(cfg(ns))?;
        }
    }
    for c in 0..3 {
        unsafe {
            s.launch_builder(&ik.sub_hinv).arg(&mut out[c]).arg(&w.tmp[c]).arg(&anu).arg(&n).arg(&nh).launch(cfg(ns))?;
        }
    }
    // force = (dx³/a) Σ g
    let lc = LaunchConfig { grid_dim: (REDUCE_BLOCKS, 1, 1), block_dim: (BLOCK, 1, 1), shared_mem_bytes: 3 * BLOCK * 8 };
    unsafe {
        s.launch_builder(&ik.sum3).arg(&g[0]).arg(&g[1]).arg(&g[2]).arg(&mut w.reduce_f64).arg(&nr64).launch(lc)?;
    }
    let parts = s.clone_dtoh(&w.reduce_f64)?;
    let mut f = [0.0f64; 3];
    for p in parts.chunks(3) {
        for c in 0..3 {
            f[c] += p[c];
        }
    }
    let dv = (2.0 * std::f64::consts::PI / n as f64).powi(3);
    dev.g_prev[stage] = g;
    w.penal = Some(pen);
    Ok((iters, rnorm / bnorm, f.map(|v| v * dv / a)))
}

impl CudaSolver {
    pub(crate) fn try_step_imex(&mut self, h: f64) -> Result<()> {
        let tab = spectral_core::imex::ars343();
        let (ns, nr) = (self.w.ns(), self.w.nr());
        if self.imex.is_none() {
            self.imex = Some(Box::new(ImexDev::new(&self.w.stream, ns, nr)?));
        }
        let CudaSolver { w, u, st, imex, ik, forcing, nu, body_force, surf, nk, .. } = self;
        let bf = *body_force;
        let dev = imex.as_mut().unwrap();
        // accumulators live in the IMEX workspace (moved out to split borrows)
        let mut racc = std::mem::replace(&mut dev.racc, [alloc3(&w.stream, 1)?, alloc3(&w.stream, 1)?, alloc3(&w.stream, 1)?]);
        let mut fin = std::mem::replace(&mut dev.fin, alloc3(&w.stream, 1)?);
        let [ka, kb, kc] = &mut racc;
        let kd = &mut fin;
        let nu = *nu;
        let f = forcing.as_ref();
        let s = w.stream.clone();

        // stage 1 (explicit): N₀ = N(uⁿ); R_m = uⁿ + h a^E_{m0} N₀; final = uⁿ
        w.rhs(u, &mut dev.nl, f, bf, None)?;
        // the physical velocity just built gives max|u| for the next CFL step
        let umax = w.max_speed2_of_up()?.sqrt();
        for (m, rm) in [(1usize, &mut *ka), (2, &mut *kb), (3, &mut *kc)] {
            for c in 0..3 {
                s.memcpy_dtod(&u[c], &mut rm[c])?;
            }
            accum(w, ik, rm, [&dev.nl, &dev.nl], [h * tab.ae[m][0], 0.0])?;
        }
        for c in 0..3 {
            s.memcpy_dtod(&u[c], &mut kd[c])?;
        }
        for i in 1..4 {
            let a = h * tab.ai[i][i];
            let ri: &[CudaSlice<float2>; 3] = match i {
                1 => ka,
                2 => kb,
                _ => kc,
            };
            let (iters, resid, force) = implicit_solve(w, ik, dev, nu, ri, a, st, i - 1, bf.is_some(), surf.as_deref_mut().map(|sd| (sd, &*nk)))?;
            dev.last_iters[i - 1] = iters;
            dev.last_resid[i - 1] = resid;
            if i == 3 {
                dev.last_force = force;
                if let Some(sd) = surf.as_deref_mut() {
                    sd.last_force = force;
                }
            }
            let inv_a = (1.0 / a) as f32;
            for c in 0..3 {
                unsafe {
                    s.launch_builder(&ik.stage_lin)
                        .arg(&mut dev.lin[c]).arg(&st[c]).arg(&ri[c]).arg(&inv_a).arg(&(ns as i64))
                        .launch(cfg(ns))?;
                }
            }
            w.rhs(st, &mut dev.nl, f, bf, None)?;
            // R_m += h(a^E_{mi} N_i + a^I_{mi} L_i) for future stages m
            for m in (i + 1)..4 {
                let rm: &mut [CudaSlice<float2>; 3] = if m == 2 { kb } else { kc };
                accum(w, ik, rm, [&dev.nl, &dev.lin], [h * tab.ae[m][i], h * tab.ai[m][i]])?;
            }
            // final += h b_i (N_i + L_i)
            accum(w, ik, kd, [&dev.nl, &dev.lin], [h * tab.b[i], h * tab.b[i]])?;
        }
        std::mem::swap(u, kd);
        dev.racc = racc;
        dev.fin = fin;
        self.last_umax = Some(umax);
        self.t += h;
        Ok(())
    }

    /// CG iterations and residuals of the last IMEX step's implicit stages.
    pub fn last_cg(&self) -> Option<([usize; 3], [f64; 3])> {
        self.imex.as_ref().map(|d| (d.last_iters, d.last_resid))
    }
}
