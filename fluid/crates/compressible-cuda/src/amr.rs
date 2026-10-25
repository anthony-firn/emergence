//! Adaptive multiresolution grids on the GPU.
//!
//! The hierarchy (which blocks exist, regridding decisions) stays with the
//! CPU `compressible::amr::Amr`; the solution lives on the device. All
//! blocks — leaves and ancestors — sit in one pool with face ghost slabs,
//! and every stage runs as a few batched kernels over the pool:
//!
//!   restriction of the edge ancestors (finest level first; fine patch
//!   gather + one 6-point pass per axis with the CPU's stencil shifts) →
//!   ghost slabs (same-level copy, 6th-order prediction from the coarser
//!   level, boundary conditions incl. characteristic) → primitives →
//!   per level and axis: sensor, midpoint Riemann fluxes, face fluxes,
//!   divergence, boundary-face store → refluxing → SSP-RK3 update.
//!
//! On a regrid step the leaves are downloaded, the CPU regrids (wavelet
//! details, refine/coarsen with 2:1 balance), and the pool is rebuilt and
//! uploaded if anything changed. In f64 the result matches the CPU
//! `Amr::step` to rounding (summation order differs in the predictions);
//! `examples/gpu_amr_validate` checks it.

use crate::{cfg, load_kernels, GpuSystem, Kern, Real, Result, RATE_BLOCKS};
use crate::Geo;
use compressible::amr::{Amr, AmrLayout};
use compressible::scheme::{FluxMode, ShockScheme};
use compressible::solver::Bc;
use compressible::System;
use cudarc::driver::{CudaSlice, CudaStream, DevicePtr, DeviceRepr, LaunchConfig, PushKernelArg};
use std::sync::Arc;

/// Device mirror of `AmrT` in kernels.cu (raw pointers into slices owned by
/// `GpuAmr`; one stream, so no cross-stream hazards).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
struct AmrT {
    blevel: u64,
    bidx: u64,
    mapoff: u64,
    mapd: u64,
    far: u64,
    bcell: i64,
    ntot: i64,
    nb0: [i32; 3],
    act: [i32; 3],
    bd: [i32; 3],
    ext: [i32; 3],
    mom: [i32; 3],
    bc: [[i32; 2]; 3],
    ng: i32,
    levels: i32,
    pad: i32,
}
unsafe impl DeviceRepr for AmrT {}

/// ancestors restricted per launch group (bounds the patch scratch)
const RCHUNK: usize = 32;

struct Pool<R: Real> {
    lay: AmrLayout,
    nblk: usize,
    /// leaves in slot order (by level); per level (first slot, count)
    leaves: Vec<usize>,
    lslots: Vec<(usize, usize)>,
    geo: Geo,
    t: AmrT,
    u: CudaSlice<R>,
    u0: CudaSlice<R>,
    w: CudaSlice<R>,
    rhs: CudaSlice<R>,
    th: CudaSlice<R>,
    mf: CudaSlice<R>,
    mu: CudaSlice<R>,
    f: CudaSlice<R>,
    uf: CudaSlice<R>,
    bl: CudaSlice<i32>,
    _blevel: CudaSlice<i32>,
    _bidx: CudaSlice<i32>,
    _mapoff: CudaSlice<i64>,
    _mapd: CudaSlice<i32>,
    edge: Vec<Option<CudaSlice<i32>>>,
    ra: CudaSlice<R>,
    rb: CudaSlice<R>,
    oka: CudaSlice<u8>,
    okb: CudaSlice<u8>,
    bf: CudaSlice<R>,
    bft: i64,
    nfc: i64,
    reflux: Vec<(usize, usize, usize, CudaSlice<i32>)>,
    /// per axis: face slabs predicted from the coarser level (8 ints per
    /// region: block, side, fine origin, coarse patch origin), count
    preg: Vec<Option<(CudaSlice<i32>, usize)>>,
    pa: CudaSlice<R>,
    pb: CudaSlice<R>,
}

pub struct GpuAmr<R: Real, const N: usize, S: System<N>> {
    /// hierarchy and host copy of the solution (current after `download`)
    pub cpu: Amr<N, S>,
    pub sys: GpuSystem,
    stream: Arc<CudaStream>,
    k: Kern,
    p: Pool<R>,
    _far: CudaSlice<R>,
    far_ptr: u64,
    wp: CudaSlice<R>,
    wts: CudaSlice<R>,
    invdx: CudaSlice<R>,
    counters: CudaSlice<u64>,
    err: CudaSlice<u64>,
    rates: CudaSlice<R>,
    stage_dt: f64,
    /// device time per phase (s, measured with synchronisation when
    /// `profile` is set): [restrict, ghosts, fluxes, reflux+update, regrid]
    pub prof: [f64; 5],
    /// regrid breakdown (s): [details, download, refine, upload, coarsen]
    pub rprof: [f64; 5],
    pub profile: bool,
}

fn ptr<T>(s: &CudaSlice<T>, st: &CudaStream) -> u64 {
    let (p, _g) = s.device_ptr(st);
    p as u64
}

fn lagrange6(x0: f64, x: f64) -> [f64; 6] {
    std::array::from_fn(|k| {
        let xk = x0 + k as f64;
        (0..6).filter(|&j| j != k).map(|j| (x - (x0 + j as f64)) / (xk - (x0 + j as f64))).product()
    })
}

impl<R: Real, const N: usize, S: System<N>> GpuAmr<R, N, S> {
    pub fn new(cpu: Amr<N, S>, sys: GpuSystem) -> Result<Self> {
        assert_eq!(N, sys.nv(), "system size mismatch");
        if cpu.relax.is_some() {
            return Err("per-cell relaxation is not available on the GPU".into());
        }
        if cpu.shock != ShockScheme::Wcns {
            return Err("only the WCNS shock path is ported to the GPU".into());
        }
        if cpu.subcycle {
            return Err("subcycling is not ported to the GPU yet".into());
        }
        let (stream, k) = load_kernels::<R>(sys)?;
        let mut farv = vec![R::default(); 6 * N];
        for d in 0..3 {
            for s in 0..2 {
                if let Bc::Dirichlet(v) | Bc::Characteristic(v) = cpu.bc[d][s] {
                    for q in 0..N {
                        farv[(d * 2 + s) * N + q] = R::from_f64(v[q]);
                    }
                }
            }
        }
        let far = stream.clone_htod(&farv)?;
        let far_ptr = ptr(&far, &stream);
        let lay = cpu.layout();
        let wpv: Vec<R> = lay.wp.iter().flat_map(|r| r.iter().map(|v| R::from_f64(*v))).collect();
        let shifts = [0i64, -1, 1, -2, 2];
        let wtsv: Vec<R> = (0..5).flat_map(|m| if m == 0 { lay.wr } else { lagrange6(-2.5 + shifts[m] as f64, 0.0) }).map(R::from_f64).collect();
        let invdxv: Vec<R> = lay.dx.iter().flat_map(|d| d.iter().map(|v| R::from_f64(1.0 / v))).collect();
        let wp = stream.clone_htod(&wpv)?;
        let wts = stream.clone_htod(&wtsv)?;
        let invdx = stream.clone_htod(&invdxv)?;
        let p = Self::build_pool(&cpu, sys, &stream, far_ptr)?;
        Ok(GpuAmr {
            counters: stream.alloc_zeros::<u64>(4)?,
            err: stream.alloc_zeros::<u64>(1)?,
            rates: stream.alloc_zeros::<R>(2 * RATE_BLOCKS as usize)?,
            cpu,
            sys,
            stream,
            k,
            p,
            _far: far,
            far_ptr,
            wp,
            wts,
            invdx,
            stage_dt: 0.0,
            prof: [0.0; 5],
            rprof: [0.0; 5],
            profile: false,
        })
    }

    fn build_pool(cpu: &Amr<N, S>, sys: GpuSystem, stream: &Arc<CudaStream>, far_ptr: u64) -> Result<Pool<R>> {
        let lay = cpu.layout();
        let act = lay.active;
        let ng = lay.ng;
        let bd = lay.bd;
        let ext: [usize; 3] = std::array::from_fn(|d| if act[d] { bd[d] + 2 * ng } else { 1 });
        let bcell = ext[0] * ext[1] * ext[2];
        let nint = bd[0] * bd[1] * bd[2];
        let nblk = lay.level.len();
        let ntot = nblk * bcell;
        // pool upload (interiors; ghosts are filled on the device)
        let mut host = vec![R::default(); N * ntot];
        for (id, b) in cpu.blocks.iter().enumerate() {
            for i in 0..bd[0] {
                for j in 0..bd[1] {
                    for kk in 0..bd[2] {
                        let c = [i, j, kk];
                        let off: usize = id * bcell + (0..3).fold(0, |o, d| o * ext[d] + if act[d] { c[d] + ng } else { 0 });
                        let v = &b.u[(i * bd[1] + j) * bd[2] + kk];
                        for q in 0..N {
                            host[q * ntot + off] = R::from_f64(v[q]);
                        }
                    }
                }
            }
        }
        let u = stream.clone_htod(&host)?;
        let u0 = stream.alloc_zeros::<R>(N * ntot)?;
        let w = stream.alloc_zeros::<R>(N * ntot)?;
        // leaves by level
        let mut leaves = Vec::new();
        let mut lslots = Vec::new();
        for l in 0..lay.levels {
            let s0 = leaves.len();
            leaves.extend((0..nblk).filter(|&i| lay.leaf[i] && lay.level[i] == l));
            lslots.push((s0, leaves.len() - s0));
        }
        let nleaf = leaves.len();
        let maxl = lslots.iter().map(|x| x.1).max().unwrap_or(1).max(1);
        let rhs = stream.alloc_zeros::<R>(N * nleaf.max(1) * nint)?;
        let nf_max = (0..3).filter(|&d| act[d]).map(|d| (bd[d] + 1) * nint / bd[d]).max().unwrap_or(1);
        let nm_max = (0..3).filter(|&d| act[d]).map(|d| (bd[d] + 5) * nint / bd[d]).max().unwrap_or(1);
        let th = stream.alloc_zeros::<R>(maxl * nf_max)?;
        let f = stream.alloc_zeros::<R>(N * maxl * nf_max)?;
        let uf = stream.alloc_zeros::<R>(maxl * nf_max)?;
        let mf = stream.alloc_zeros::<R>(N * maxl * nm_max)?;
        let mu = stream.alloc_zeros::<R>(maxl * nm_max)?;
        let bl = stream.clone_htod(&leaves.iter().map(|&i| i as i32).collect::<Vec<_>>())?;
        let blevel = stream.clone_htod(&lay.level.iter().map(|&l| l as i32).collect::<Vec<_>>())?;
        let bidx = stream.clone_htod(&lay.idx.iter().flat_map(|c| c.iter().map(|&v| v as i32)).collect::<Vec<_>>())?;
        let mut offs = Vec::new();
        let mut md = Vec::new();
        for lv in &lay.map {
            offs.push(md.len() as i64);
            md.extend_from_slice(lv);
        }
        let mapoff = stream.clone_htod(&offs)?;
        let mapd = stream.clone_htod(&md)?;
        let edge: Vec<Option<CudaSlice<i32>>> =
            lay.edge.iter().map(|e| if e.is_empty() { Ok(None) } else { stream.clone_htod(&e.iter().map(|&i| i as i32).collect::<Vec<_>>()).map(Some) }).collect::<std::result::Result<_, _>>()?;
        let psz: usize = (0..3).map(|d| if act[d] { 2 * bd[d] + 8 } else { 1 }).product();
        let maxe = lay.edge.iter().map(|e| e.len()).max().unwrap_or(0).min(RCHUNK).max(1);
        let ra = stream.alloc_zeros::<R>(N * maxe * psz)?;
        let rb = stream.alloc_zeros::<R>(N * maxe * psz)?;
        let oka = stream.alloc_zeros::<u8>(maxe * psz)?;
        let okb = stream.alloc_zeros::<u8>(maxe * psz)?;
        let nfc = (0..3).filter(|&d| act[d]).map(|d| nint / bd[d]).max().unwrap_or(1);
        let bft = nleaf.max(1) * 6 * nfc;
        let bf = stream.alloc_zeros::<R>(N * bft)?;
        // refluxing tasks per (axis, side), slots for block ids
        let mut slot = vec![usize::MAX; nblk];
        for (s, &id) in leaves.iter().enumerate() {
            slot[id] = s;
        }
        let mut reflux = Vec::new();
        for d in 0..3 {
            for side in 0..2 {
                let mut tk: Vec<i32> = Vec::new();
                for rf in lay.reflux.iter().filter(|r| r.d == d && r.side == side) {
                    let mut row = [0i32; 13];
                    row[0] = slot[rf.block] as i32;
                    row[1] = rf.cell as i32;
                    row[2] = rf.face as i32;
                    row[3] = lay.level[rf.block] as i32;
                    row[4] = rf.fine.len() as i32;
                    for (j, &(fb, ff)) in rf.fine.iter().enumerate() {
                        row[5 + 2 * j] = slot[fb] as i32;
                        row[6 + 2 * j] = ff as i32;
                    }
                    tk.extend_from_slice(&row);
                }
                if !tk.is_empty() {
                    let n = tk.len() / 13;
                    reflux.push((d, side, n, stream.clone_htod(&tk)?));
                }
            }
        }
        // face slabs of leaves whose face neighbour region is coarser
        assert!((0..3).all(|d| !act[d] || bd[d] % 2 == 0), "block size must be even");
        // coarse patch along the slab axis: ⌊ng/2⌋ + 7 cells, all used for odd ng
        assert!(ng % 2 == 1, "ghost depth must be odd");
        let mut preg = Vec::new();
        let mut pmax = 1usize;
        for d in 0..3 {
            let mut rows: Vec<i32> = Vec::new();
            if act[d] {
                for &id in &leaves {
                    let l = lay.level[id];
                    if l == 0 {
                        continue;
                    }
                    let idx = lay.idx[id];
                    for side in 0..2usize {
                        let nbk = lay.nb0[d] << l;
                        let mut nk = idx;
                        nk[d] += if side == 1 { 1 } else { -1 };
                        if nk[d] < 0 || nk[d] >= nbk {
                            if cpu.bc[d][side] != Bc::Periodic {
                                continue;
                            }
                            nk[d] = nk[d].rem_euclid(nbk);
                        }
                        let nbl: [i64; 3] = std::array::from_fn(|e| if act[e] { lay.nb0[e] << l } else { 1 });
                        let slot = (nk[0] * nbl[1] + nk[1]) * nbl[2] + nk[2];
                        if lay.map[l][slot as usize] >= 0 {
                            continue;
                        }
                        let mut row = [0i32; 8];
                        row[0] = id as i32;
                        row[1] = side as i32;
                        for e in 0..3 {
                            if !act[e] {
                                continue;
                            }
                            let f0 = if e == d {
                                if side == 1 { (idx[e] + 1) * bd[e] as i64 } else { idx[e] * bd[e] as i64 - ng as i64 }
                            } else {
                                idx[e] * bd[e] as i64
                            };
                            row[2 + e] = f0 as i32;
                            row[5 + e] = (f0.div_euclid(2) - 3) as i32;
                        }
                        rows.extend_from_slice(&row);
                    }
                }
            }
            let n = rows.len() / 8;
            if n > 0 {
                let stage: usize = (0..3).map(|e| if !act[e] { 1 } else if e == d { ng / 2 + 7 } else { (bd[e] / 2 + 6).max(bd[e]) }).product();
                pmax = pmax.max(n * stage);
                preg.push(Some((stream.clone_htod(&rows)?, n)));
            } else {
                preg.push(None);
            }
        }
        let pa = stream.alloc_zeros::<R>(N * pmax)?;
        let pb = stream.alloc_zeros::<R>(N * pmax)?;
        let code = |b: &Bc<N>| -> i32 {
            match b {
                Bc::Periodic => 0,
                Bc::Outflow => 1,
                Bc::Reflect => 2,
                Bc::Dirichlet(_) => 3,
                Bc::Characteristic(_) => 4,
            }
        };
        let t = AmrT {
            blevel: ptr(&blevel, stream),
            bidx: ptr(&bidx, stream),
            mapoff: ptr(&mapoff, stream),
            mapd: ptr(&mapd, stream),
            far: far_ptr,
            bcell: bcell as i64,
            ntot: ntot as i64,
            nb0: lay.nb0.map(|v| v as i32),
            act: act.map(|a| a as i32),
            bd: bd.map(|v| v as i32),
            ext: ext.map(|v| v as i32),
            mom: [sys.mom(0), sys.mom(1), sys.mom(2)],
            bc: std::array::from_fn(|d| [code(&cpu.bc[d][0]), code(&cpu.bc[d][1])]),
            ng: ng as i32,
            levels: lay.levels as i32,
            pad: 0,
        };
        let geo = Geo {
            m: ext.map(|v| v as i32),
            n: bd.map(|v| v as i32),
            g: std::array::from_fn(|d| if act[d] { ng as i32 } else { 0 }),
            _pad: 0,
            stride: [(ext[1] * ext[2]) as i64, ext[2] as i64, 1],
            ncell: bcell as i64,
            ntot: ntot as i64,
            rtot: (nleaf.max(1) * nint) as i64,
            nb: 0,
            s0: 0,
        };
        Ok(Pool {
            lay,
            nblk,
            leaves,
            lslots,
            geo,
            t,
            u,
            u0,
            w,
            rhs,
            th,
            mf,
            mu,
            f,
            uf,
            bl,
            _blevel: blevel,
            _bidx: bidx,
            _mapoff: mapoff,
            _mapd: mapd,
            edge,
            ra,
            rb,
            oka,
            okb,
            bf,
            bft: bft as i64,
            nfc: nfc as i64,
            reflux,
            preg,
            pa,
            pb,
        })
    }

    fn eos(&self) -> ([R; 8], i32) {
        let (a, kap) = self.sys.eos_args();
        (a.map(R::from_f64), kap)
    }

    fn act(&self) -> [bool; 3] {
        self.p.lay.active
    }

    fn tick(&self) -> Result<std::time::Instant> {
        if self.profile {
            self.stream.synchronize()?;
        }
        Ok(std::time::Instant::now())
    }
    fn tock(&mut self, t0: std::time::Instant, i: usize) -> Result<()> {
        if self.profile {
            self.stream.synchronize()?;
            self.prof[i] += t0.elapsed().as_secs_f64();
        }
        Ok(())
    }

    /// Restrict the edge ancestors, finest level first.
    fn restrict(&mut self) -> Result<()> {
        let (e, kap) = self.eos();
        let act = self.act();
        let bd = self.p.lay.bd;
        let psz: [usize; 3] = std::array::from_fn(|d| if act[d] { 2 * bd[d] + 8 } else { 1 });
        let dims: Vec<usize> = (0..3).filter(|&d| act[d]).collect();
        let t = self.p.t;
        for l in (0..self.p.lay.levels.saturating_sub(1)).rev() {
            let Some(list) = self.p.edge[l].as_ref() else { continue };
            let total = self.p.lay.edge[l].len();
            let mut c0 = 0;
            while c0 < total {
                let nt = (total - c0).min(RCHUNK);
                let tl = list.slice(c0..c0 + nt);
                let (nti, li) = (nt as i32, l as i32);
                {
                    let mut lb = self.stream.launch_builder(&self.k.amr_rgather);
                    for v in e.iter() {
                        lb.arg(v);
                    }
                    lb.arg(&kap).arg(&t).arg(&self.p.u).arg(&tl).arg(&nti).arg(&li).arg(&mut self.p.ra).arg(&mut self.p.oka);
                    unsafe { lb.launch(cfg(nt * psz.iter().product::<usize>()))? };
                }
                let mut sz = psz;
                for (pi, &d) in dims.iter().enumerate() {
                    let fin = (pi + 1 == dims.len()) as i32;
                    let mut so = sz;
                    so[d] = bd[d];
                    let (s0, s1, s2, di) = (sz[0] as i32, sz[1] as i32, sz[2] as i32, d as i32);
                    // ping-pong: even passes ra → rb, odd passes rb → ra
                    let mut lb = self.stream.launch_builder(&self.k.amr_rpass);
                    lb.arg(&t).arg(&tl).arg(&nti).arg(&di).arg(&s0).arg(&s1).arg(&s2);
                    if pi % 2 == 0 {
                        lb.arg(&self.p.ra).arg(&self.p.oka).arg(&mut self.p.rb).arg(&mut self.p.okb);
                    } else {
                        lb.arg(&self.p.rb).arg(&self.p.okb).arg(&mut self.p.ra).arg(&mut self.p.oka);
                    }
                    lb.arg(&self.wts).arg(&fin).arg(&mut self.p.u);
                    unsafe { lb.launch(cfg(nt * so.iter().product::<usize>()))? };
                    sz = so;
                }
                c0 += nt;
            }
        }
        Ok(())
    }

    fn ghosts(&mut self) -> Result<()> {
        let (e, kap) = self.eos();
        let act = self.act();
        let bd = self.p.lay.bd;
        let ng = self.p.lay.ng;
        let nb = self.p.leaves.len();
        let t = self.p.t;
        for d in 0..3 {
            if !act[d] {
                continue;
            }
            let per = 2 * ng * bd[0] * bd[1] * bd[2] / bd[d];
            let (nbi, di) = (nb as i32, d as i32);
            let mut lb = self.stream.launch_builder(&self.k.amr_ghost);
            for v in e.iter() {
                lb.arg(v);
            }
            let sep: i32 = 1;
            lb.arg(&kap).arg(&t).arg(&mut self.p.u).arg(&self.p.bl).arg(&nbi).arg(&di).arg(&self.wp).arg(&mut self.err).arg(&sep);
            unsafe { lb.launch(cfg(per * nb))? };
        }
        // coarse–fine slabs: gather coarse patches, one pass per axis
        let dims: Vec<usize> = (0..3).filter(|&e| act[e]).collect();
        for d in 0..3 {
            let Some((reg, nreg)) = self.p.preg[d].as_ref() else { continue };
            let nreg = *nreg;
            let szc: [usize; 3] = std::array::from_fn(|e| if !act[e] { 1 } else if e == d { ng / 2 + 7 } else { bd[e] / 2 + 6 });
            let fine: [usize; 3] = std::array::from_fn(|e| if !act[e] { 1 } else if e == d { ng } else { bd[e] });
            let nri = nreg as i32;
            {
                let (c0, c1, c2) = (szc[0] as i32, szc[1] as i32, szc[2] as i32);
                let mut lb = self.stream.launch_builder(&self.k.amr_pgather);
                for v in e.iter() {
                    lb.arg(v);
                }
                lb.arg(&kap).arg(&t).arg(&self.p.u).arg(reg).arg(&nri).arg(&c0).arg(&c1).arg(&c2).arg(&mut self.p.pa).arg(&mut self.err);
                unsafe { lb.launch(cfg(nreg * szc.iter().product::<usize>()))? };
            }
            let mut sz = szc;
            for (pi, &a) in dims.iter().enumerate() {
                let fin = (pi + 1 == dims.len()) as i32;
                let mut so = sz;
                so[a] = fine[a];
                let (i0, i1, i2, ai, oa) = (sz[0] as i32, sz[1] as i32, sz[2] as i32, a as i32, fine[a] as i32);
                let mut lb = self.stream.launch_builder(&self.k.amr_ppass);
                lb.arg(&t).arg(reg).arg(&nri).arg(&ai).arg(&i0).arg(&i1).arg(&i2).arg(&oa);
                if pi % 2 == 0 {
                    lb.arg(&self.p.pa).arg(&mut self.p.pb);
                } else {
                    lb.arg(&self.p.pb).arg(&mut self.p.pa);
                }
                lb.arg(&self.wp).arg(&fin).arg(&mut self.p.u);
                unsafe { lb.launch(cfg(nreg * so.iter().product::<usize>()))? };
                sz = so;
            }
        }
        Ok(())
    }

    fn prim(&mut self) -> Result<()> {
        let (e, kap) = self.eos();
        let mut geo = self.p.geo;
        geo.nb = self.p.leaves.len() as i32;
        let mut lb = self.stream.launch_builder(&self.k.prim);
        for v in e.iter() {
            lb.arg(v);
        }
        lb.arg(&kap).arg(&self.p.u).arg(&mut self.p.w).arg(&geo).arg(&self.p.bl);
        unsafe { lb.launch(cfg(self.p.leaves.len() * geo.ncell as usize))? };
        Ok(())
    }

    /// Ghost-filled primitives of all leaves (restriction first).
    fn prepare(&mut self) -> Result<()> {
        let t0 = self.tick()?;
        self.restrict()?;
        self.tock(t0, 0)?;
        let t0 = self.tick()?;
        self.ghosts()?;
        self.prim()?;
        self.tock(t0, 1)
    }

    fn rhs(&mut self) -> Result<()> {
        self.prepare()?;
        let t0 = self.tick()?;
        let (e, kap) = self.eos();
        let act = self.act();
        let ndim = act.iter().filter(|&&a| a).count() as f64;
        let mode: i32 = match self.cpu.mode {
            FluxMode::Hybrid => 0,
            FluxMode::Weno => 1,
            FluxMode::Central => 2,
        };
        let tol = R::from_f64(self.cpu.tol);
        let bd = self.p.lay.bd;
        let nint = bd[0] * bd[1] * bd[2];
        for l in 0..self.p.lay.levels {
            let (s0, cnt) = self.p.lslots[l];
            if cnt == 0 {
                continue;
            }
            let mut geo = self.p.geo;
            geo.nb = cnt as i32;
            geo.s0 = s0 as i32;
            let bl = self.p.bl.slice(s0..s0 + cnt);
            let mut first: i32 = 1;
            for d in 0..3 {
                if !act[d] {
                    continue;
                }
                let dx = self.p.lay.dx[l][d];
                let lambda = R::from_f64(if self.cpu.positivity { ndim * self.stage_dt / dx } else { 0.0 });
                let di = d as i32;
                let nf = (bd[d] + 1) * nint / bd[d];
                let nm = (bd[d] + 5) * nint / bd[d];
                {
                    let mut lb = self.stream.launch_builder(&self.k.theta);
                    for v in e.iter() {
                        lb.arg(v);
                    }
                    lb.arg(&kap).arg(&self.p.w).arg(&geo).arg(&bl).arg(&di).arg(&mode).arg(&tol).arg(&mut self.p.th);
                    unsafe { lb.launch(cfg(nf * cnt))? };
                }
                if mode != 2 {
                    let mut lb = self.stream.launch_builder(&self.k.midflux);
                    for v in e.iter() {
                        lb.arg(v);
                    }
                    lb.arg(&kap).arg(&self.p.w).arg(&geo).arg(&bl).arg(&di).arg(&self.p.th).arg(&mut self.p.mf).arg(&mut self.p.mu).arg(&mut self.counters);
                    unsafe { lb.launch(cfg(nm * cnt))? };
                }
                {
                    let mut lb = self.stream.launch_builder(&self.k.flux);
                    for v in e.iter() {
                        lb.arg(v);
                    }
                    lb.arg(&kap)
                        .arg(&self.p.w)
                        .arg(&geo)
                        .arg(&bl)
                        .arg(&di)
                        .arg(&lambda)
                        .arg(&self.p.th)
                        .arg(&self.p.mf)
                        .arg(&self.p.mu)
                        .arg(&mut self.p.f)
                        .arg(&mut self.p.uf)
                        .arg(&mut self.counters);
                    unsafe { lb.launch(cfg(nf * cnt))? };
                }
                {
                    let dxr = R::from_f64(dx);
                    let mut lb = self.stream.launch_builder(&self.k.accum);
                    for v in e.iter() {
                        lb.arg(v);
                    }
                    lb.arg(&kap).arg(&self.p.w).arg(&geo).arg(&bl).arg(&di).arg(&dxr).arg(&self.p.f).arg(&self.p.uf).arg(&mut self.p.rhs).arg(&first);
                    unsafe { lb.launch(cfg(nint * cnt))? };
                }
                if !self.p.reflux.is_empty() {
                    let mut lb = self.stream.launch_builder(&self.k.amr_faces);
                    lb.arg(&geo).arg(&di).arg(&self.p.f).arg(&mut self.p.bf).arg(&self.p.bft).arg(&self.p.nfc);
                    unsafe { lb.launch(cfg(2 * nint / bd[d] * cnt))? };
                }
                first = 0;
            }
        }
        self.tock(t0, 2)?;
        let t0 = self.tick()?;
        let rtot = self.p.geo.rtot;
        let nint_l = nint as i64;
        for (d, side, n, tk) in &self.p.reflux {
            let (ni, di, si) = (*n as i32, *d as i32, *side as i32);
            let mut lb = self.stream.launch_builder(&self.k.amr_reflux);
            lb.arg(&ni).arg(tk).arg(&self.invdx).arg(&di).arg(&si).arg(&self.p.bf).arg(&self.p.bft).arg(&self.p.nfc).arg(&mut self.p.rhs).arg(&rtot).arg(&nint_l);
            unsafe { lb.launch(cfg(*n))? };
        }
        self.tock(t0, 3)
    }

    fn rk(&mut self, al: f64, be: f64, dt: f64) -> Result<()> {
        let mut geo = self.p.geo;
        geo.nb = self.p.leaves.len() as i32;
        geo.s0 = 0;
        let bd = self.p.lay.bd;
        let (a, b, t) = (R::from_f64(al), R::from_f64(be), R::from_f64(dt));
        let mut lb = self.stream.launch_builder(&self.k.rk);
        lb.arg(&mut self.p.u).arg(&self.p.u0).arg(&self.p.rhs).arg(&geo).arg(&self.p.bl).arg(&a).arg(&b).arg(&t);
        unsafe { lb.launch(cfg(bd[0] * bd[1] * bd[2] * self.p.leaves.len()))? };
        Ok(())
    }

    /// One global SSP-RK3 step; regrids on the CPU when due.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        self.stage_dt = dt;
        self.stream.memcpy_dtod(&self.p.u, &mut self.p.u0)?;
        self.rhs()?;
        self.rk(0.0, 1.0, dt)?;
        self.rhs()?;
        self.rk(0.75, 0.25, dt)?;
        self.rhs()?;
        self.rk(1.0 / 3.0, 2.0 / 3.0, dt)?;
        self.cpu.advance_clock(dt);
        if self.cpu.regrid_is_due() {
            let t0 = std::time::Instant::now();
            let changed = self.regrid()?;
            self.cpu.schedule_regrid(changed);
            self.prof[4] += t0.elapsed().as_secs_f64();
        }
        Ok(())
    }

    /// Wavelet details of all leaves on the device (block id order, 0 for
    /// ancestors), from freshly restricted ancestors and filled ghosts.
    fn details(&mut self) -> Result<Vec<f64>> {
        self.prepare()?;
        let (e, kap) = self.eos();
        let nl = self.p.leaves.len();
        let mut det = self.stream.alloc_zeros::<R>(nl.max(1))?;
        let wr = self.stream.clone_htod(&self.p.lay.wr.map(R::from_f64))?;
        let geo = self.p.geo;
        let mut lb = self.stream.launch_builder(&self.k.amr_detail);
        for v in e.iter() {
            lb.arg(v);
        }
        lb.arg(&kap).arg(&geo).arg(&self.p.w).arg(&self.p.bl).arg(&wr).arg(&mut det);
        unsafe { lb.launch(LaunchConfig { grid_dim: (nl as u32, 1, 1), block_dim: (256, 1, 1), shared_mem_bytes: 0 })? };
        let h = self.stream.clone_dtoh(&det)?;
        let mut out = vec![0.0; self.p.nblk];
        for (s, &id) in self.p.leaves.iter().enumerate() {
            out[id] = h[s].to_f64();
        }
        Ok(out)
    }

    /// `Amr::regrid` with the details computed on the device; the CPU only
    /// edits the tree (refined children are predicted, coarsened parents
    /// restricted on the host from the downloaded solution).
    fn regrid(&mut self) -> Result<bool> {
        let clock = std::time::Instant::now;
        let before = self.cpu.structure();
        let t0 = clock();
        let det = self.details()?; // also restricts the ancestors
        self.stream.synchronize()?;
        self.rprof[0] += t0.elapsed().as_secs_f64();
        let t0 = clock();
        self.download()?;
        self.rprof[1] += t0.elapsed().as_secs_f64();
        let t0 = clock();
        let flags = self.cpu.flags_from_details(&det);
        let protected = self.cpu.regrid_refine(&flags);
        let refined = self.cpu.structure() != before;
        self.rprof[2] += t0.elapsed().as_secs_f64();
        let flags = if refined {
            let t0 = clock();
            self.p = Self::build_pool(&self.cpu, self.sys, &self.stream, self.far_ptr)?;
            self.stream.synchronize()?;
            self.rprof[3] += t0.elapsed().as_secs_f64();
            let t0 = clock();
            let det = self.details()?;
            self.rprof[0] += t0.elapsed().as_secs_f64();
            self.cpu.flags_from_details(&det)
        } else {
            flags
        };
        let mid = self.cpu.structure();
        let t0 = clock();
        self.cpu.regrid_coarsen(&flags, &protected);
        self.rprof[4] += t0.elapsed().as_secs_f64();
        if self.cpu.structure() != mid {
            let t0 = clock();
            self.p = Self::build_pool(&self.cpu, self.sys, &self.stream, self.far_ptr)?;
            self.stream.synchronize()?;
            self.rprof[3] += t0.elapsed().as_secs_f64();
        }
        Ok(self.cpu.structure() != before)
    }

    /// CFL (and multiphase dilatation) time step over the leaves.
    pub fn dt_cfl(&mut self) -> Result<f64> {
        self.prepare()?;
        let (e, kap) = self.eos();
        let (mut r1, mut r2) = (0.0f64, 0.0f64);
        for l in 0..self.p.lay.levels {
            let (s0, cnt) = self.p.lslots[l];
            if cnt == 0 {
                continue;
            }
            let mut geo = self.p.geo;
            geo.nb = cnt as i32;
            geo.s0 = s0 as i32;
            let bl = self.p.bl.slice(s0..s0 + cnt);
            let dx = self.p.lay.dx[l].map(R::from_f64);
            let lc = LaunchConfig { grid_dim: (RATE_BLOCKS, 1, 1), block_dim: (256, 1, 1), shared_mem_bytes: 0 };
            let mut lb = self.stream.launch_builder(&self.k.rate);
            for v in e.iter() {
                lb.arg(v);
            }
            lb.arg(&kap).arg(&self.p.w).arg(&geo).arg(&bl).arg(&dx[0]).arg(&dx[1]).arg(&dx[2]).arg(&mut self.rates);
            unsafe { lb.launch(lc)? };
            let r = self.stream.clone_dtoh(&self.rates)?;
            for b in 0..RATE_BLOCKS as usize {
                r1 = r1.max(r[2 * b].to_f64());
                r2 = r2.max(r[2 * b + 1].to_f64());
            }
        }
        let err = self.stream.clone_dtoh(&self.err)?[0];
        if err > 0 {
            return Err(format!("{err} ghost cells could not be filled (hierarchy not 2:1 balanced?)").into());
        }
        let mut dt = self.cpu.cfl / r1;
        if r2 > 0.0 {
            dt = dt.min(0.5 / r2);
        }
        Ok(dt)
    }

    pub fn run_to(&mut self, t_end: f64) -> Result<usize> {
        let mut n = 0;
        while self.cpu.t < t_end * (1.0 - 1e-14) {
            let dt = self.dt_cfl()?.min(t_end - self.cpu.t);
            if !(dt.is_finite() && dt > 0.0) {
                return Err("time step collapsed".into());
            }
            self.step(dt)?;
            n += 1;
        }
        Ok(n)
    }

    /// Copy the device solution into `cpu.blocks` (all blocks).
    pub fn download(&mut self) -> Result<()> {
        let h = self.stream.clone_dtoh(&self.p.u)?;
        let lay = &self.p.lay;
        let (act, bd, ng) = (lay.active, lay.bd, lay.ng);
        let ext: [usize; 3] = std::array::from_fn(|d| if act[d] { bd[d] + 2 * ng } else { 1 });
        let bcell = ext[0] * ext[1] * ext[2];
        let ntot = self.p.nblk * bcell;
        for (id, b) in self.cpu.blocks.iter_mut().enumerate() {
            for i in 0..bd[0] {
                for j in 0..bd[1] {
                    for kk in 0..bd[2] {
                        let c = [i, j, kk];
                        let off: usize = id * bcell + (0..3).fold(0, |o, d| o * ext[d] + if act[d] { c[d] + ng } else { 0 });
                        let v = &mut b.u[(i * bd[1] + j) * bd[2] + kk];
                        for q in 0..N {
                            v[q] = h[q * ntot + off].to_f64();
                        }
                    }
                }
            }
        }
        Ok(())
    }

    pub fn t(&self) -> f64 {
        self.cpu.t
    }

    /// [shock-path faces, first-order fallbacks, positivity-limited faces]
    pub fn counters(&self) -> Result<[u64; 3]> {
        let c = self.stream.clone_dtoh(&self.counters)?;
        Ok([c[0], c[1], c[2]])
    }

    pub fn leaf_count(&self) -> usize {
        self.p.leaves.len()
    }
}
