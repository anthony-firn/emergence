//! GPU port of the a posteriori error certificate (`spectral_core::certify`).
//! Same algorithm: exact truncation residual of v × ω on a padded grid,
//! Hermite time residual at the step midpoint, μ = max(−λ_min(sym ∇v)) on
//! the padded grid, and the weak–strong bound y' ≤ 2μy + ‖R‖²_{H⁻¹}/(2ν).
//! Reductions are in f64; fields are f32, which puts a floor of roughly
//! 1e-7·‖u‖ under the residuals (an additional floating-point caveat).

use super::*;
use spectral_core::certify::CertStep;
use spectral_core::grid::wavenumber;

struct CertKernels {
    trunc: CudaFunction,
    strain: CudaFunction,
    min_eig: CudaFunction,
    sub_visc: CudaFunction,
    hermite: CudaFunction,
    spec_sub: CudaFunction,
    hm1: CudaFunction,
    scatter: CudaFunction,
    scale: CudaFunction,
}

pub struct CudaCertifier {
    k: CertKernels,
    m: i32,
    mh: i32,
    kmax: i32,
    r2c_p: CudaFft,
    c2r_p: CudaFft,
    pn: CudaSlice<i32>,
    pm: CudaSlice<i32>,
    ones: CudaSlice<f32>,
    npairs: i32,
    pspec: [CudaSlice<float2>; 3],
    pw: [CudaSlice<float2>; 3],
    ptmp: CudaSlice<float2>,
    pphys: [CudaSlice<f32>; 6],
    prev_u: [CudaSlice<float2>; 3],
    prev_f: [CudaSlice<float2>; 3],
    f1: [CudaSlice<float2>; 3],
    vm: [CudaSlice<float2>; 3],
    vdot: [CudaSlice<float2>; 3],
    fm: [CudaSlice<float2>; 3],
    started: bool,
    prev_mu: f64,
    prev_r2: f64,
    /// certified bound on ‖u − v‖² (box-mean L²)
    pub y: f64,
    pub history: Vec<CertStep>,
}

impl CudaCertifier {
    /// `initial_error`: a bound on ‖u(0) − v(0)‖ (0 when the initial data is
    /// exactly representable). Requires a 2/3-rule grid.
    pub fn new(solver: &CudaSolver, initial_error: f64) -> Result<Self> {
        let g = &solver.grid;
        if g.spherical {
            return Err("the certificate assumes the exact 2/3-rule discretization".into());
        }
        let kmax = g.kmax();
        let mut m = (3 * g.n) / 2;
        if m < 4 * kmax + 2 {
            m = 4 * kmax + 2;
        }
        m += m % 2;
        let mh = m / 2 + 1;
        let (mut pn, mut pmv) = (Vec::new(), Vec::new());
        for ni in 0..g.spec_len() {
            if !g.kept(ni) {
                continue;
            }
            let (i, j, l) = g.spec_ijl(ni);
            let w = |k: i64| k.rem_euclid(m as i64) as usize;
            pn.push(ni as i32);
            pmv.push(((w(wavenumber(i, g.n)) * m + w(wavenumber(j, g.n))) * mh + l) as i32);
        }
        let s = &solver.w.stream;
        let module = &solver.module;
        let k = CertKernels {
            trunc: module.load_function("trunc_resid_partial")?,
            strain: module.load_function("strain_spec")?,
            min_eig: module.load_function("min_eig_partial")?,
            sub_visc: module.load_function("sub_viscous")?,
            hermite: module.load_function("hermite_mid")?,
            spec_sub: module.load_function("spec_sub")?,
            hm1: module.load_function("hm1_partial")?,
            scatter: module.load_function("spec_scatter")?,
            scale: module.load_function("spec_scale")?,
        };
        let mi = m as i32;
        let (nsp, nrp, ns) = (m * m * mh, m * m * m, g.spec_len());
        Ok(Self {
            k,
            m: mi,
            mh: mh as i32,
            kmax: kmax as i32,
            r2c_p: CudaFft::plan_3d(mi, mi, mi, cufftType::CUFFT_R2C, s.clone())?,
            c2r_p: CudaFft::plan_3d(mi, mi, mi, cufftType::CUFFT_C2R, s.clone())?,
            npairs: pn.len() as i32,
            ones: s.clone_htod(&vec![1.0f32; pn.len()])?,
            pn: s.clone_htod(&pn)?,
            pm: s.clone_htod(&pmv)?,
            pspec: alloc3(s, nsp)?,
            pw: alloc3(s, nsp)?,
            ptmp: s.alloc_zeros(nsp)?,
            pphys: [s.alloc_zeros(nrp)?, s.alloc_zeros(nrp)?, s.alloc_zeros(nrp)?, s.alloc_zeros(nrp)?, s.alloc_zeros(nrp)?, s.alloc_zeros(nrp)?],
            prev_u: alloc3(s, ns)?,
            prev_f: alloc3(s, ns)?,
            f1: alloc3(s, ns)?,
            vm: alloc3(s, ns)?,
            vdot: alloc3(s, ns)?,
            fm: alloc3(s, ns)?,
            started: false,
            prev_mu: 0.0,
            prev_r2: 0.0,
            y: initial_error * initial_error,
            history: Vec::new(),
        })
    }

    fn reduce_sum(s: &Arc<CudaStream>, buf: &CudaSlice<f64>) -> Result<f64> {
        Ok(s.clone_dtoh(buf)?[..REDUCE_BLOCKS as usize].iter().sum())
    }

    /// (‖P(v×ω)|discarded‖²_{H⁻¹}, max −λ_min(sym ∇v)) for N-grid spectrum `uh`.
    fn evaluate(&mut self, w: &mut Work, uh: &[CudaSlice<float2>; 3]) -> Result<(f64, f64)> {
        let s = w.stream.clone();
        let (m, mh) = (self.m, self.mh);
        let nsp = (m as usize) * (m as usize) * (mh as usize);
        let nrp = (m as usize).pow(3);
        let nr64 = nrp as i64;
        let zf = 0.0f32;
        let pc = LaunchConfig { grid_dim: ((self.npairs as u32).div_ceil(256), 1, 1), block_dim: (256, 1, 1), shared_mem_bytes: 0 };
        for c in 0..3 {
            s.memset_zeros(&mut self.pspec[c])?;
            unsafe {
                s.launch_builder(&self.k.scatter)
                    .arg(&mut self.pspec[c]).arg(&uh[c]).arg(&self.pn).arg(&self.pm).arg(&self.ones).arg(&self.npairs)
                    .launch(pc)?;
            }
            s.memcpy_dtod(&self.pspec[c], &mut self.ptmp)?;
            self.c2r_p.exec_c2r(&mut self.ptmp, &mut self.pphys[c])?;
        }
        {
            let [w0, w1, w2] = &mut self.pw;
            unsafe {
                s.launch_builder(&w.k.curl)
                    .arg(&self.pspec[0]).arg(&self.pspec[1]).arg(&self.pspec[2])
                    .arg(&mut *w0).arg(&mut *w1).arg(&mut *w2)
                    .arg(&m).arg(&mh).arg(&zf).arg(&zf).arg(&zf)
                    .launch(cfg(nsp))?;
            }
        }
        for c in 0..3 {
            self.c2r_p.exec_c2r(&mut self.pw[c], &mut self.pphys[3 + c])?;
        }
        {
            let [u0, u1, u2, w0, w1, w2] = &mut self.pphys;
            unsafe {
                s.launch_builder(&w.k.cross)
                    .arg(&*u0).arg(&*u1).arg(&*u2).arg(&mut *w0).arg(&mut *w1).arg(&mut *w2).arg(&nr64)
                    .launch(cfg(nrp))?;
            }
        }
        let inv = 1.0f32 / nrp as f32;
        for c in 0..3 {
            self.r2c_p.exec_r2c(&self.pphys[3 + c], &mut self.pw[c])?;
            unsafe {
                s.launch_builder(&self.k.scale).arg(&mut self.pw[c]).arg(&inv).arg(&(nsp as i64)).launch(cfg(nsp))?;
            }
        }
        let lc = |bytes: u32| LaunchConfig { grid_dim: (REDUCE_BLOCKS, 1, 1), block_dim: (BLOCK, 1, 1), shared_mem_bytes: bytes };
        unsafe {
            s.launch_builder(&self.k.trunc)
                .arg(&self.pw[0]).arg(&self.pw[1]).arg(&self.pw[2])
                .arg(&mut w.reduce_f64).arg(&m).arg(&mh).arg(&self.kmax)
                .launch(lc(BLOCK * 8))?;
        }
        let r2 = Self::reduce_sum(&s, &w.reduce_f64)?;
        // strain components S_ij on the padded grid
        let comps = [(0i32, 0i32), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];
        for (q, &(i, j)) in comps.iter().enumerate() {
            unsafe {
                s.launch_builder(&self.k.strain)
                    .arg(&self.pspec[i as usize]).arg(&self.pspec[j as usize]).arg(&mut self.ptmp)
                    .arg(&i).arg(&j).arg(&m).arg(&mh)
                    .launch(cfg(nsp))?;
            }
            self.c2r_p.exec_c2r(&mut self.ptmp, &mut self.pphys[q])?;
        }
        unsafe {
            let p = &self.pphys;
            s.launch_builder(&self.k.min_eig)
                .arg(&p[0]).arg(&p[1]).arg(&p[2]).arg(&p[3]).arg(&p[4]).arg(&p[5])
                .arg(&mut w.reduce_f64).arg(&nr64)
                .launch(lc(BLOCK * 8))?;
        }
        let parts = s.clone_dtoh(&w.reduce_f64)?;
        let mu = parts[..REDUCE_BLOCKS as usize].iter().cloned().fold(0.0f64, f64::max);
        Ok((r2, mu))
    }
}

impl CudaSolver {
    /// Full tendency F(v) = P[(v × ω)^] − ν|k|² v̂ into `out`.
    fn tendency(w: &mut Work, k: &CertKernels, nu: f64, uh: &[CudaSlice<float2>; 3], out: &mut [CudaSlice<float2>; 3]) -> Result<()> {
        w.rhs(uh, out, None, None, None)?;
        let (n, nh, nuf) = (w.n, w.nh, nu as f32);
        for c in 0..3 {
            unsafe {
                w.stream.launch_builder(&k.sub_visc).arg(&mut out[c]).arg(&uh[c]).arg(&nuf).arg(&n).arg(&nh).launch(cfg(w.ns()))?;
            }
        }
        Ok(())
    }

    /// Advance one IF-RK4 step and extend the certificate (unforced,
    /// unpenalized, 2/3-rule runs only).
    pub fn step_certified(&mut self, h: f64, cert: &mut CudaCertifier) -> Result<CertStep> {
        if self.forcing.is_some() || self.w.penal.is_some() || self.surf.is_some() || self.body_force.is_some() || self.spherical {
            return Err("the certificate covers unforced 2/3-rule Navier–Stokes only".into());
        }
        let nu = self.nu;
        let s = self.w.stream.clone();
        let ns = self.w.ns();
        if !cert.started {
            let (r2, mu) = cert.evaluate(&mut self.w, &self.u)?;
            Self::tendency(&mut self.w, &cert.k, nu, &self.u, &mut cert.prev_f)?;
            for c in 0..3 {
                s.memcpy_dtod(&self.u[c], &mut cert.prev_u[c])?;
            }
            cert.prev_mu = mu;
            cert.prev_r2 = r2;
            cert.started = true;
        }
        self.try_step(h)?;
        Self::tendency(&mut self.w, &cert.k, nu, &self.u, &mut cert.f1)?;
        let hf = h as f32;
        for c in 0..3 {
            unsafe {
                s.launch_builder(&cert.k.hermite)
                    .arg(&mut cert.vm[c]).arg(&mut cert.vdot[c])
                    .arg(&cert.prev_u[c]).arg(&self.u[c]).arg(&cert.prev_f[c]).arg(&cert.f1[c])
                    .arg(&hf).arg(&(ns as i64))
                    .launch(cfg(ns))?;
            }
        }
        Self::tendency(&mut self.w, &cert.k, nu, &cert.vm, &mut cert.fm)?;
        // R_time = P(v̇ − F(v)) on retained modes, in H⁻¹
        for c in 0..3 {
            unsafe {
                s.launch_builder(&cert.k.spec_sub).arg(&mut cert.vdot[c]).arg(&cert.fm[c]).arg(&(ns as i64)).launch(cfg(ns))?;
            }
        }
        let (n, nh) = (self.w.n, self.w.nh);
        let (zero_i, one, zf, kc2) = (0i32, 1.0f32, 0.0f32, 0.0f32);
        {
            let [v0, v1, v2] = &mut cert.vdot;
            unsafe {
                s.launch_builder(&self.w.k.project)
                    .arg(&mut *v0).arg(&mut *v1).arg(&mut *v2)
                    .arg(&self.w.tmp[0]).arg(&self.w.tmp[1]).arg(&self.w.tmp[2])
                    .arg(&zero_i).arg(&one).arg(&n).arg(&nh)
                    .arg(&zero_i).arg(&zf).arg(&zf).arg(&zf)
                    .arg(&kc2).arg(&zf).arg(&zf).arg(&zf)
                    .launch(cfg(ns))?;
            }
        }
        let lc = LaunchConfig { grid_dim: (REDUCE_BLOCKS, 1, 1), block_dim: (BLOCK, 1, 1), shared_mem_bytes: BLOCK * 8 };
        unsafe {
            s.launch_builder(&cert.k.hm1)
                .arg(&cert.vdot[0]).arg(&cert.vdot[1]).arg(&cert.vdot[2])
                .arg(&mut self.w.reduce_f64).arg(&n).arg(&nh)
                .launch(lc)?;
        }
        let rtime2 = CudaCertifier::reduce_sum(&s, &self.w.reduce_f64)?;
        let vm = std::mem::replace(&mut cert.vm, alloc3(&s, 1)?);
        let (rtr_m, mu_m) = cert.evaluate(&mut self.w, &vm)?;
        cert.vm = vm;
        let (rtr_1, mu_1) = cert.evaluate(&mut self.w, &self.u)?;

        let r2_mid = rtime2 + rtr_m;
        let int_r2 = h / 6.0 * (cert.prev_r2 + 4.0 * r2_mid + rtr_1);
        let mu = cert.prev_mu.max(mu_m).max(mu_1);
        cert.y = (2.0 * mu * h).exp() * (cert.y + int_r2 / (2.0 * nu));
        for c in 0..3 {
            s.memcpy_dtod(&self.u[c], &mut cert.prev_u[c])?;
        }
        std::mem::swap(&mut cert.prev_f, &mut cert.f1);
        cert.prev_mu = mu_1;
        cert.prev_r2 = rtr_1;
        let st = CertStep { t: self.t, resid_trunc: rtr_m.sqrt(), resid_time: rtime2.sqrt(), mu, bound: cert.y.sqrt() };
        cert.history.push(st);
        Ok(st)
    }
}
