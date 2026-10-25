//! GPU surface-constrained walls (point forces), mirroring
//! `spectral_core::imex::surface_solve`: Gaussian-gridding NUFFT on the
//! device, Lagrange-multiplier CG on the host (markers are few).

use super::*;
use spectral_core::grid::wavenumber;
use spectral_core::Surface;

/// Gaussian half width on the GPU: ~2e-6 relative accuracy, matched to f32.
const HALF_WIDTH: usize = 6;

pub(crate) struct NufftKernels {
    scatter: CudaFunction,
    gather: CudaFunction,
    interp: CudaFunction,
    spread: CudaFunction,
    pub(crate) add_hinv: CudaFunction,
    ext_build: CudaFunction,
    ext_adjoint: CudaFunction,
    spec_scale: CudaFunction,
}

impl NufftKernels {
    pub(crate) fn load(m: &Arc<cudarc::driver::CudaModule>) -> Result<Self> {
        Ok(Self {
            scatter: m.load_function("spec_scatter")?,
            gather: m.load_function("spec_gather")?,
            interp: m.load_function("nufft_interp")?,
            spread: m.load_function("nufft_spread")?,
            add_hinv: m.load_function("add_hinv_scaled")?,
            ext_build: m.load_function("ext_build")?,
            ext_adjoint: m.load_function("ext_adjoint")?,
            spec_scale: m.load_function("spec_scale")?,
        })
    }
}

pub(crate) struct DevSurface {
    m: i32,
    tau: f32,
    h: f32,
    r2c_m: CudaFft,
    c2r_m: CudaFft,
    pn: CudaSlice<i32>,
    pm: CudaSlice<i32>,
    inv_phi: CudaSlice<f32>,
    npairs: i32,
    grid_r: CudaSlice<f32>,
    grid_c: CudaSlice<float2>,
    pts: CudaSlice<f32>,
    pub(crate) nm: usize,
    pub(crate) u_s: [f64; 3],
    /// per-marker prescribed velocities (overrides u_s), as in `Surface::target`
    pub(crate) targets: Option<Vec<[f64; 3]>>,
    vals: CudaSlice<f32>,
    forces: CudaSlice<f32>,
    /// warm starts per implicit stage, 3·nm each
    pub(crate) f_prev: Vec<Vec<f64>>,
    pub(crate) last_force: [f64; 3],
    /// smooth forcing extension (least squares), if requested
    pub(crate) ext: Option<DevExtension>,
    /// host copy (part ranges for per-part forces)
    pub(crate) surface: Surface,
}

/// Device copy of the sparse extension operator (CSR by shell grid point).
pub(crate) struct DevExtension {
    grid_idx: CudaSlice<i32>,
    row_ptr: CudaSlice<i32>,
    col: CudaSlice<i32>,
    val: CudaSlice<f32>,
    nrows: i32,
    nb: usize,
    colsum: Vec<f64>,
    amp: CudaSlice<f32>,
    res: CudaSlice<f32>,
    a_prev: Vec<Vec<f64>>,
}

impl DevSurface {
    /// Force on each part (last implicit stage).
    pub(crate) fn part_forces(&self) -> Vec<[f64; 3]> {
        match &self.ext {
            Some(e) => {
                let op = spectral_core::surface::ExtensionOp { rows: Vec::new(), colsum: e.colsum.clone() };
                self.surface.part_forces_extension(&op, &e.a_prev[2])
            }
            None => self.surface.part_forces_point(&self.f_prev[2]),
        }
    }

    pub(crate) fn set_targets(&mut self, t: Vec<[f64; 3]>) {
        assert_eq!(t.len(), self.nm);
        self.targets = Some(t);
    }

    #[inline]
    fn target(&self, m: usize, c: usize) -> f64 {
        match &self.targets {
            Some(t) => t[m][c],
            None => self.u_s[c],
        }
    }

    pub(crate) fn new(stream: &Arc<CudaStream>, grid: &Grid, surface: &Surface) -> Result<Self> {
        let kmax = grid.kmax();
        let modes = 2 * kmax + 1;
        let mut m = 2 * modes;
        m += m % 2;
        let r = m as f64 / modes as f64;
        let hw = HALF_WIDTH as f64;
        let tau = std::f64::consts::PI * hw / ((1.0 - 1.0 / r).sqrt() * (m * m) as f64);
        let phi = |k: f64| 2.0 * (std::f64::consts::PI * tau).sqrt() * (-k * k * tau).exp();
        let mh = m / 2 + 1;
        let (mut pn, mut pmv, mut ip) = (Vec::new(), Vec::new(), Vec::new());
        for ni in 0..grid.spec_len() {
            if !grid.kept(ni) {
                continue;
            }
            let (i, j, l) = grid.spec_ijl(ni);
            let (kx, ky, kz) = (wavenumber(i, grid.n), wavenumber(j, grid.n), l as i64);
            let w = |k: i64| k.rem_euclid(m as i64) as usize;
            pn.push(ni as i32);
            pmv.push(((w(kx) * m + w(ky)) * mh + l) as i32);
            ip.push((1.0 / (phi(kx as f64) * phi(ky as f64) * phi(kz as f64))) as f32);
        }
        let mi = m as i32;
        let pts: Vec<f32> = surface.points.iter().flat_map(|p| p.iter().map(|&v| v as f32)).collect();
        let nm = surface.points.len();
        Ok(Self {
            m: mi,
            tau: tau as f32,
            h: (2.0 * std::f64::consts::PI / m as f64) as f32,
            r2c_m: CudaFft::plan_3d(mi, mi, mi, cufftType::CUFFT_R2C, stream.clone())?,
            c2r_m: CudaFft::plan_3d(mi, mi, mi, cufftType::CUFFT_C2R, stream.clone())?,
            npairs: pn.len() as i32,
            pn: stream.clone_htod(&pn)?,
            pm: stream.clone_htod(&pmv)?,
            inv_phi: stream.clone_htod(&ip)?,
            grid_r: stream.alloc_zeros(m * m * m)?,
            grid_c: stream.alloc_zeros(m * m * mh)?,
            pts: stream.clone_htod(&pts)?,
            nm,
            u_s: surface.u_solid,
            targets: surface.u_markers.clone(),
            vals: stream.alloc_zeros(nm)?,
            forces: stream.alloc_zeros(nm)?,
            f_prev: vec![vec![0.0; 3 * nm]; 3],
            last_force: [0.0; 3],
            surface: surface.clone(),
            ext: match surface.extension {
                None => None,
                Some(e) => {
                    let op = surface.extension_op(grid, e, surface.spacing);
                    let nb = surface.centres().len();
                    let (mut gi, mut rp, mut col, mut val) = (Vec::new(), vec![0i32], Vec::new(), Vec::new());
                    for (pidx, coup) in &op.rows {
                        gi.push(*pidx as i32);
                        for &(m, w) in coup {
                            col.push(m as i32);
                            val.push(w as f32);
                        }
                        rp.push(col.len() as i32);
                    }
                    Some(DevExtension {
                        nrows: gi.len() as i32,
                        grid_idx: stream.clone_htod(&gi)?,
                        row_ptr: stream.clone_htod(&rp)?,
                        col: stream.clone_htod(&col)?,
                        val: stream.clone_htod(&val)?,
                        nb,
                        colsum: op.colsum.clone(),
                        amp: stream.alloc_zeros(nb)?,
                        res: stream.alloc_zeros(nb)?,
                        a_prev: vec![vec![0.0; 3 * nb]; 3],
                    })
                }
            },
        })
    }

    fn pcfg(&self, n: usize) -> LaunchConfig {
        LaunchConfig { grid_dim: ((n as u32).div_ceil(128), 1, 1), block_dim: (128, 1, 1), shared_mem_bytes: 0 }
    }

    /// Values at the markers of the field with N-grid coefficients `src`.
    pub(crate) fn interpolate(&mut self, s: &Arc<CudaStream>, k: &NufftKernels, src: &CudaSlice<float2>) -> Result<Vec<f64>> {
        s.memset_zeros(&mut self.grid_c)?;
        let lc = self.pcfg(self.npairs as usize);
        unsafe {
            s.launch_builder(&k.scatter)
                .arg(&mut self.grid_c).arg(src).arg(&self.pn).arg(&self.pm).arg(&self.inv_phi).arg(&self.npairs)
                .launch(lc)?;
        }
        self.c2r_m.exec_c2r(&mut self.grid_c, &mut self.grid_r)?;
        let (nmi, hw) = (self.nm as i32, HALF_WIDTH as i32);
        let lc = self.pcfg(self.nm);
        unsafe {
            s.launch_builder(&k.interp)
                .arg(&mut self.vals).arg(&self.grid_r).arg(&self.pts)
                .arg(&nmi).arg(&self.m).arg(&hw).arg(&self.tau).arg(&self.h)
                .launch(lc)?;
        }
        Ok(s.clone_dtoh(&self.vals)?.into_iter().map(f64::from).collect())
    }

    /// Retained-mode coefficients of Σ F_m δ(x − X_m) into `out` (N grid).
    pub(crate) fn spread(&mut self, s: &Arc<CudaStream>, k: &NufftKernels, forces: &[f64], out: &mut CudaSlice<float2>) -> Result<()> {
        let f32s: Vec<f32> = forces.iter().map(|&v| v as f32).collect();
        s.memcpy_htod(&f32s, &mut self.forces)?;
        s.memset_zeros(&mut self.grid_r)?;
        let (nmi, hw) = (self.nm as i32, HALF_WIDTH as i32);
        let lc = self.pcfg(self.nm);
        unsafe {
            s.launch_builder(&k.spread)
                .arg(&mut self.grid_r).arg(&self.pts).arg(&self.forces)
                .arg(&nmi).arg(&self.m).arg(&hw).arg(&self.tau).arg(&self.h)
                .launch(lc)?;
        }
        self.r2c_m.exec_r2c(&self.grid_r, &mut self.grid_c)?;
        s.memset_zeros(out)?;
        let extra = 1.0f32 / (self.m as f32).powi(3);
        let lc = self.pcfg(self.npairs as usize);
        unsafe {
            s.launch_builder(&k.gather)
                .arg(out).arg(&self.grid_c).arg(&self.pn).arg(&self.pm).arg(&self.inv_phi).arg(&extra).arg(&self.npairs)
                .launch(lc)?;
        }
        Ok(())
    }
}

/// Surface constraint in a GPU implicit stage. On entry `out` = H⁻¹r̂;
/// scratch spectra `sa`, `sb` are free IMEX buffers. Mirrors the CPU solve.
#[allow(clippy::too_many_arguments)]
pub(crate) fn surface_solve(
    w: &mut Work,
    nk: &NufftKernels,
    sd: &mut DevSurface,
    nu: f64,
    a: f64,
    out: &mut [CudaSlice<float2>; 3],
    sa: &mut [CudaSlice<float2>; 3],
    sb: &mut [CudaSlice<float2>; 3],
    stage: usize,
    keep_mean: bool,
) -> Result<(usize, f64, [f64; 3])> {
    let s = w.stream.clone();
    let nm = sd.nm;
    let (n, nh) = (w.n, w.nh);
    let ns = w.ns();
    let (anu, af) = ((a * nu) as f32, a as f32);
    let keep: i32 = keep_mean as i32;
    let (zero_i, one, zf) = (0i32, 1.0f32, 0.0f32);
    let mut b = vec![0.0; 3 * nm];
    let mut wall_scale2 = 0.0;
    for c in 0..3 {
        let v = sd.interpolate(&s, nk, &out[c])?;
        for m in 0..nm {
            let tg = sd.target(m, c);
            b[c * nm + m] = tg - v[m];
            wall_scale2 += v[m] * v[m] + tg * tg;
        }
    }
    // spread F into `sa`, project, and leave a·H⁻¹·P·S F in `sb`
    let apply_into = |f: &[f64], sa: &mut [CudaSlice<float2>; 3], sb: &mut [CudaSlice<float2>; 3], sd: &mut DevSurface| -> Result<()> {
        for c in 0..3 {
            sd.spread(&s, nk, &f[c * nm..(c + 1) * nm], &mut sa[c])?;
        }
        {
            let [s0, s1, s2] = sa;
            unsafe {
                s.launch_builder(&w.k.project)
                    .arg(&mut *s0).arg(&mut *s1).arg(&mut *s2)
                    .arg(&w.tmp[0]).arg(&w.tmp[1]).arg(&w.tmp[2])
                    .arg(&zero_i).arg(&one).arg(&n).arg(&nh)
                    .arg(&keep).arg(&zf).arg(&zf).arg(&zf)
                    .arg(&w.kc2).arg(&zf).arg(&zf).arg(&zf)
                    .launch(cfg(ns))?;
            }
        }
        for c in 0..3 {
            s.memset_zeros(&mut sb[c])?;
            unsafe {
                s.launch_builder(&nk.add_hinv).arg(&mut sb[c]).arg(&sa[c]).arg(&af).arg(&anu).arg(&n).arg(&nh).launch(cfg(ns))?;
            }
        }
        Ok(())
    };
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p * q).sum::<f64>();
    let mut f = std::mem::take(&mut sd.f_prev[stage]);
    let mut ap = vec![0.0; 3 * nm];
    let matvec = |x: &[f64], ax: &mut [f64], sa: &mut [CudaSlice<float2>; 3], sb: &mut [CudaSlice<float2>; 3], sd: &mut DevSurface| -> Result<()> {
        apply_into(x, sa, sb, sd)?;
        for c in 0..3 {
            let v = sd.interpolate(&s, nk, &sb[c])?;
            ax[c * nm..(c + 1) * nm].copy_from_slice(&v);
        }
        Ok(())
    };
    matvec(&f, &mut ap, sa, sb, sd)?;
    let mut r: Vec<f64> = b.iter().zip(&ap).map(|(p, q)| p - q).collect();
    let mut p = r.clone();
    // Converge relative to the wall-velocity scale, not ‖b‖: with small time
    // steps b is tiny and an f32 solve cannot reach a tolerance relative to it.
    let bnorm = wall_scale2.sqrt().max(dot(&b, &b).sqrt()).max(1e-30);
    let mut rr = dot(&r, &r);
    let mut iters = 0;
    while rr.sqrt() > CG_TOL_SURF * bnorm && iters < 300 {
        matvec(&p, &mut ap, sa, sb, sd)?;
        let alpha = rr / dot(&p, &ap);
        f.iter_mut().zip(&p).for_each(|(x, pv)| *x += alpha * pv);
        r.iter_mut().zip(&ap).for_each(|(x, av)| *x -= alpha * av);
        let rr_new = dot(&r, &r);
        let beta = rr_new / rr;
        rr = rr_new;
        p.iter_mut().zip(&r).for_each(|(pv, rv)| *pv = rv + beta * *pv);
        iters += 1;
    }
    // Û = H⁻¹ r̂ + a H⁻¹ P S F
    apply_into(&f, sa, sb, sd)?;
    for c in 0..3 {
        // out += sb (add_hinv with zero viscosity is a plain axpy)
        unsafe {
            s.launch_builder(&nk.add_hinv).arg(&mut out[c]).arg(&sb[c]).arg(&one).arg(&zf).arg(&n).arg(&nh).launch(cfg(ns))?;
        }
    }
    let force = [0, 1, 2].map(|c| -f[c * nm..(c + 1) * nm].iter().sum::<f64>());
    sd.f_prev[stage] = f;
    Ok((iters, rr.sqrt() / bnorm, force))
}

/// CG tolerance for the GPU surface solve (f32 NUFFT ~2e-6).
const CG_TOL_SURF: f64 = 1e-5;

/// Smooth forcing extension in a GPU implicit stage (least squares, CGLS on
/// the host). Mirrors `spectral_core::imex::extension_solve`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn extension_solve(
    w: &mut Work,
    nk: &NufftKernels,
    sd: &mut DevSurface,
    nu: f64,
    a: f64,
    out: &mut [CudaSlice<float2>; 3],
    sa: &mut [CudaSlice<float2>; 3],
    sb: &mut [CudaSlice<float2>; 3],
    stage: usize,
    keep_mean: bool,
) -> Result<(usize, f64, [f64; 3])> {
    let s = w.stream.clone();
    let nm = sd.nm;
    let mut ext = sd.ext.take().expect("extension");
    let nb = ext.nb;
    let (n, nh) = (w.n, w.nh);
    let (ns, nr) = (w.ns(), w.nr());
    let (anu, af) = ((a * nu) as f32, a as f32);
    let keep: i32 = keep_mean as i32;
    let (zero_i, one, zf) = (0i32, 1.0f32, 0.0f32);
    let inv_n3 = 1.0f32 / nr as f32;
    let two_pi3 = (2.0 * std::f64::consts::PI).powi(3) as f32;
    let rows_cfg = LaunchConfig { grid_dim: ((ext.nrows as u32).div_ceil(128), 1, 1), block_dim: (128, 1, 1), shared_mem_bytes: 0 };

    let project = |w: &Work, v: &mut [CudaSlice<float2>; 3], scale: f32| -> Result<()> {
        let [v0, v1, v2] = v;
        unsafe {
            s.launch_builder(&w.k.project)
                .arg(&mut *v0).arg(&mut *v1).arg(&mut *v2)
                .arg(&w.tmp[0]).arg(&w.tmp[1]).arg(&w.tmp[2])
                .arg(&zero_i).arg(&scale).arg(&n).arg(&nh)
                .arg(&keep).arg(&zf).arg(&zf).arg(&zf)
                .arg(&w.kc2).arg(&zf).arg(&zf).arg(&zf)
                .launch(cfg(ns))?;
        }
        Ok(())
    };
    // E·A into w.wp (physical), forward to `sa`, projected: sa = P F(E A).
    // Returns the box means of E·A (read before projection, which zeroes
    // k = 0 when the mean flow is held fixed) only when `want_mean`.
    let build_project = |amp: &[f64], w: &mut Work, sa: &mut [CudaSlice<float2>; 3], ext: &mut DevExtension, want_mean: bool| -> Result<[f64; 3]> {
        let mut means = [0.0; 3];
        for c in 0..3 {
            let a32: Vec<f32> = amp[c * nb..(c + 1) * nb].iter().map(|&v| v as f32).collect();
            s.memcpy_htod(&a32, &mut ext.amp)?;
            s.memset_zeros(&mut w.wp[c])?;
            unsafe {
                s.launch_builder(&nk.ext_build)
                    .arg(&mut w.wp[c]).arg(&ext.amp).arg(&ext.grid_idx).arg(&ext.row_ptr).arg(&ext.col).arg(&ext.val).arg(&ext.nrows)
                    .launch(rows_cfg)?;
            }
            w.r2c.exec_r2c(&w.wp[c], &mut sa[c])?;
            if want_mean {
                means[c] = s.clone_dtoh(&sa[c].slice(0..1))?[0].x as f64 * inv_n3 as f64;
            }
        }
        project(w, sa, inv_n3)?;
        Ok(means)
    };
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).map(|(p, q)| p * q).sum::<f64>();
    let mut b = vec![0.0; 3 * nm];
    let mut wall_scale2 = 0.0;
    for c in 0..3 {
        let v = sd.interpolate(&s, nk, &out[c])?;
        for m in 0..nm {
            let tg = sd.target(m, c);
            b[c * nm + m] = tg - v[m];
            wall_scale2 += v[m] * v[m] + tg * tg;
        }
    }
    let fwd = |amp: &[f64], ax: &mut [f64], w: &mut Work, sa: &mut [CudaSlice<float2>; 3], sb: &mut [CudaSlice<float2>; 3], sd: &mut DevSurface, ext: &mut DevExtension| -> Result<()> {
        build_project(amp, w, sa, ext, false)?;
        for c in 0..3 {
            s.memset_zeros(&mut sb[c])?;
            unsafe {
                s.launch_builder(&nk.add_hinv).arg(&mut sb[c]).arg(&sa[c]).arg(&af).arg(&anu).arg(&n).arg(&nh).launch(cfg(ns))?;
            }
            let v = sd.interpolate(&s, nk, &sb[c])?;
            ax[c * nm..(c + 1) * nm].copy_from_slice(&v);
        }
        Ok(())
    };
    let adj = |y: &[f64], res: &mut [f64], w: &mut Work, sa: &mut [CudaSlice<float2>; 3], sb: &mut [CudaSlice<float2>; 3], sd: &mut DevSurface, ext: &mut DevExtension| -> Result<()> {
        for c in 0..3 {
            sd.spread(&s, nk, &y[c * nm..(c + 1) * nm], &mut sa[c])?;
            unsafe {
                s.launch_builder(&nk.spec_scale).arg(&mut sa[c]).arg(&two_pi3).arg(&(ns as i64)).launch(cfg(ns))?;
            }
        }
        project(w, sa, one)?;
        let coef = af * inv_n3;
        for c in 0..3 {
            s.memset_zeros(&mut sb[c])?;
            unsafe {
                s.launch_builder(&nk.add_hinv).arg(&mut sb[c]).arg(&sa[c]).arg(&coef).arg(&anu).arg(&n).arg(&nh).launch(cfg(ns))?;
            }
            w.c2r.exec_c2r(&mut sb[c], &mut w.wp[c])?;
            s.memset_zeros(&mut ext.res)?;
            unsafe {
                s.launch_builder(&nk.ext_adjoint)
                    .arg(&mut ext.res).arg(&w.wp[c]).arg(&ext.grid_idx).arg(&ext.row_ptr).arg(&ext.col).arg(&ext.val).arg(&ext.nrows)
                    .launch(rows_cfg)?;
            }
            let r = s.clone_dtoh(&ext.res)?;
            for (dst, v) in res[c * nb..(c + 1) * nb].iter_mut().zip(r) {
                *dst = v as f64;
            }
        }
        Ok(())
    };

    // CGLS (as on the CPU)
    let mut x = std::mem::take(&mut ext.a_prev[stage]);
    // tolerance relative to the wall-velocity scale (see surface_solve)
    let bnorm = wall_scale2.sqrt().max(dot(&b, &b).sqrt()).max(1e-30);
    let mut ax = vec![0.0; 3 * nm];
    fwd(&x, &mut ax, w, sa, sb, sd, &mut ext)?;
    let mut r: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
    let mut sv = vec![0.0; 3 * nb];
    adj(&r, &mut sv, w, sa, sb, sd, &mut ext)?;
    let mut p = sv.clone();
    let mut gamma = dot(&sv, &sv);
    let gamma0 = gamma.max(1e-300);
    let mut products = 2;
    let mut q = vec![0.0; 3 * nm];
    let mut history: Vec<f64> = Vec::new();
    while products + 2 <= EXT_MAX_PRODUCTS {
        let rnorm = dot(&r, &r).sqrt() / bnorm;
        // Stagnation: < 2% progress over 5 iterations. Warm starts usually
        // land at the reachable residual floor, and a 20-iteration window
        // spent 40 products per solve just confirming that.
        let stagnant = history.len() >= 5 && rnorm > 0.98 * history[history.len() - 5];
        if rnorm < EXT_TOL || gamma < EXT_TOL * EXT_TOL * gamma0 || stagnant {
            break;
        }
        history.push(rnorm);
        fwd(&p, &mut q, w, sa, sb, sd, &mut ext)?;
        let alpha = gamma / dot(&q, &q).max(1e-300);
        x.iter_mut().zip(&p).for_each(|(xv, pv)| *xv += alpha * pv);
        r.iter_mut().zip(&q).for_each(|(rv, qv)| *rv -= alpha * qv);
        adj(&r, &mut sv, w, sa, sb, sd, &mut ext)?;
        products += 2;
        let gamma_new = dot(&sv, &sv);
        let beta = gamma_new / gamma;
        gamma = gamma_new;
        p.iter_mut().zip(&sv).for_each(|(pv, s2)| *pv = s2 + beta * *pv);
    }
    // Û = H⁻¹r̂ + a H⁻¹ P E A ; force = −∫ E A = −V · mean(E A)
    let means = build_project(&x, w, sa, &mut ext, true)?;
    let vol = (2.0 * std::f64::consts::PI).powi(3);
    let force = means.map(|m| -m * vol);
    for c in 0..3 {
        unsafe {
            s.launch_builder(&nk.add_hinv).arg(&mut out[c]).arg(&sa[c]).arg(&af).arg(&anu).arg(&n).arg(&nh).launch(cfg(ns))?;
        }
    }
    ext.a_prev[stage] = x;
    let resid = dot(&r, &r).sqrt() / bnorm;
    sd.ext = Some(ext);
    Ok((products, resid, force))
}

const EXT_TOL: f64 = 1e-6;
const EXT_MAX_PRODUCTS: usize = 400;
