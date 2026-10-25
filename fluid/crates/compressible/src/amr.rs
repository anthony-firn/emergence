//! Adaptive multiresolution grids with a controlled error.
//!
//! The domain is covered by a block-structured octree: every leaf block is
//! a B^D uniform patch at some level l (cell size Δx₀/2^l) that runs the
//! unchanged hybrid scheme of `scheme.rs`. Which level a region needs is
//! decided by interpolating-wavelet details,
//!   d_i = u_i − P(R(u))_i ,
//! the difference between the field and its 6th-order prediction from the
//! next coarser level (R: 6th-order midpoint restriction, P: 6th-order
//! prediction). d is precisely the error made by representing the region
//! one level coarser, so thresholding at ε bounds the representation error
//! by ε; for a stable scheme the adaptive solution then stays within
//! O(n ε) of the uniform finest-grid solution after n steps (Cohen, Kaber,
//! Müller & Postel 2003). The bound is checked, not assumed:
//! `examples/amr_validate` compares against the uniform finest grid.
//!
//! * Details are measured on the system's sensor values (ln ρ, ln(p + p∞),
//!   α …), i.e. relative errors. A block refines when its largest detail
//!   exceeds ε; its 2^D children merge back only when all of theirs are
//!   below ε/128 (details of smooth data fall by 2⁶ = 64 per level for
//!   these 6th-order operators, so a smaller hysteresis would make blocks
//!   refine and coarsen alternately); one block layer of buffer around refined blocks and
//!   regridding every few steps keep moving features inside fine blocks.
//! * Neighbouring leaves differ by at most one level (2:1 balance).
//! * Ancestor blocks are kept and hold the restriction of their children
//!   (6th order where the stencil lies inside the refined region, the
//!   conservative child average at its edge), refreshed before every
//!   right-hand side. Ghost cells are filled by sampling: a block at that
//!   level (leaf or ancestor) → copy, else 6th-order prediction from the
//!   level below — the recursion only ever goes to coarser levels, so it
//!   terminates. Physical boundaries map through the boundary conditions.
//! * Conservation across levels: block faces coincide between levels, so
//!   the coarse face flux is replaced by the mean of the fine face fluxes
//!   covering it (refluxing; one global time step, no subcycling yet).

use crate::scheme::{line_fluxes, FluxMode, SchemeOpts, ShockScheme, Stats};
use crate::solver::Bc;
use crate::system::{System, V};
use rayon::prelude::*;

/// 6-point Lagrange weights at x for nodes x₀ … x₀+5.
fn lagrange6(x0: f64, x: f64) -> [f64; 6] {
    std::array::from_fn(|k| {
        let xk = x0 + k as f64;
        (0..6).filter(|&j| j != k).map(|j| (x - (x0 + j as f64)) / (xk - (x0 + j as f64))).product()
    })
}

#[derive(Clone, Copy, Debug)]
pub struct AmrConfig {
    /// cells per block side (active dimensions)
    pub b: usize,
    /// number of levels (0 … levels−1)
    pub levels: usize,
    /// detail threshold
    pub eps: f64,
    /// regrid every this many steps
    pub regrid_every: usize,
    /// blocks per side at level 0 (0 marks an inactive axis)
    pub base_blocks: [usize; 3],
    /// ghost layers (≥ 5 for the WCNS path)
    pub ghosts: usize,
}

type Key = (u8, [i64; 3]);

/// Dense per-level block index (block coordinates → block id): O(1)
/// lookups in the ghost-filling hot path instead of hashing.
#[derive(Clone, Default)]
struct BlockIndex {
    /// blocks per side at level 0
    nb0: [i64; 3],
    active: [bool; 3],
    levels: Vec<Vec<i32>>,
}

impl BlockIndex {
    fn new(nb0: [i64; 3], active: [bool; 3], levels: usize) -> Self {
        let mut b = BlockIndex { nb0, active, levels: Vec::new() };
        for l in 0..levels {
            let n: usize = (0..3).map(|d| b.nb(l, d) as usize).product();
            b.levels.push(vec![-1; n]);
        }
        b
    }
    #[inline]
    fn nb(&self, l: usize, d: usize) -> i64 {
        if self.active[d] {
            self.nb0[d] << l
        } else {
            1
        }
    }
    #[inline]
    fn slot(&self, k: &Key) -> Option<usize> {
        let l = k.0 as usize;
        if l >= self.levels.len() {
            return None;
        }
        let mut s = 0usize;
        for d in 0..3 {
            let n = self.nb(l, d);
            let c = k.1[d];
            if c < 0 || c >= n {
                return None;
            }
            s = s * n as usize + c as usize;
        }
        Some(s)
    }
    #[inline]
    fn get(&self, k: &Key) -> Option<&usize> {
        // ids are stored as i32; return a reference into a thread-local slot
        let s = self.slot(k)?;
        let v = self.levels[k.0 as usize][s];
        if v < 0 {
            None
        } else {
            Some(&IDS[v as usize])
        }
    }
    #[inline]
    fn contains_key(&self, k: &Key) -> bool {
        self.slot(k).is_some_and(|s| self.levels[k.0 as usize][s] >= 0)
    }
    fn insert(&mut self, k: Key, id: usize) {
        if let Some(s) = self.slot(&k) {
            self.levels[k.0 as usize][s] = id as i32;
        }
    }
    fn remove(&mut self, k: &Key) {
        if let Some(s) = self.slot(k) {
            self.levels[k.0 as usize][s] = -1;
        }
    }
}

/// Identity table so `BlockIndex::get` can hand out `&usize` like a map.
static IDS: [usize; 1 << 20] = {
    let mut a = [0usize; 1 << 20];
    let mut i = 0;
    while i < (1 << 20) {
        a[i] = i;
        i += 1;
    }
    a
};

impl std::ops::Index<&Key> for BlockIndex {
    type Output = usize;
    fn index(&self, k: &Key) -> &usize {
        self.get(k).expect("block not present")
    }
}

/// See `Amr::layout`.
#[derive(Clone, Debug)]
pub struct AmrLayout {
    pub ng: usize,
    pub bd: [usize; 3],
    pub active: [bool; 3],
    pub levels: usize,
    /// level-0 blocks per side
    pub nb0: [i64; 3],
    /// per block id: level, block coordinates, leaf
    pub level: Vec<usize>,
    pub idx: Vec<[i64; 3]>,
    pub leaf: Vec<bool>,
    /// per level: dense block index (row-major block coordinates → id, −1)
    pub map: Vec<Vec<i32>>,
    /// per level: ancestors read by ghost filling (border a leaf)
    pub edge: Vec<Vec<usize>>,
    /// coarse faces corrected by refluxing, in the CPU order
    pub reflux: Vec<RefluxFace>,
    /// prediction weights by parity, centred restriction weights
    pub wp: [[f64; 6]; 2],
    pub wr: [f64; 6],
    pub dx: Vec<[f64; 3]>,
}

/// One coarse boundary face whose flux is replaced by the mean of the fine
/// fluxes over it: leaf `block`, axis `d`, side, face index (row-major over
/// the two other axes), interior cell index next to it, fine (leaf, face).
#[derive(Clone, Debug)]
pub struct RefluxFace {
    pub block: usize,
    pub d: usize,
    pub side: usize,
    pub face: usize,
    pub cell: usize,
    pub fine: Vec<(usize, usize)>,
}

#[derive(Clone)]
pub struct Block<const N: usize> {
    pub level: usize,
    pub idx: [i64; 3],
    /// leaf (evolved) or ancestor (holds the restriction of its children)
    pub leaf: bool,
    pub u: Vec<V<N>>,
    u0: Vec<V<N>>,
    u_save: Vec<V<N>>,
    rhs: Vec<V<N>>,
    /// boundary face fluxes per direction and side (B^{D−1} faces each)
    face: [[Vec<V<N>>; 2]; 3],
}

pub struct Amr<const N: usize, S: System<N>> {
    pub sys: S,
    pub cfg: AmrConfig,
    pub lo: [f64; 3],
    pub hi: [f64; 3],
    pub bc: [[Bc<N>; 2]; 3],
    pub mode: FluxMode,
    pub shock: ShockScheme,
    pub tol: f64,
    pub cfl: f64,
    pub positivity: bool,
    pub t: f64,
    pub steps: usize,
    pub stats: Stats,
    pub blocks: Vec<Block<N>>,
    map: BlockIndex,
    ng: usize,
    wp: [[f64; 6]; 2],
    wr: [f64; 6],
    stage_dt: f64,
    /// wall time per phase: [restrict, ghosts, fluxes, reflux, regrid, update]
    pub prof: [f64; 6],
    /// current regrid interval (doubles after a regrid that changed nothing,
    /// up to 2× the configured interval — moving shocks must not outrun the
    /// one-block buffer; resets when the grid changes)
    regrid_interval: usize,
    next_regrid: usize,
    /// Subcycling (local time stepping): level l advances with Δt/2^l,
    /// coarse ghosts interpolated in time, flux registers keep exact
    /// conservation. `false`: one global time step for all levels.
    pub subcycle: bool,
    level_t0: Vec<f64>,
    level_dt: Vec<f64>,
    /// coarse-side and fine-side flux registers on coarse–fine faces,
    /// keyed by (coarse leaf block, direction, side, face index)
    creg: std::collections::HashMap<(usize, u8, u8, usize), V<N>>,
    freg: std::collections::HashMap<(usize, u8, u8, usize), V<N>>,
    /// per-cell relaxation after each step (as `Solver::relax`)
    pub relax: Option<Box<dyn Fn(&S, &mut V<N>, f64) -> bool + Send + Sync>>,
}

impl<const N: usize, S: System<N>> Amr<N, S> {
    fn active(&self, d: usize) -> bool {
        self.cfg.base_blocks[d] > 0 && (self.hi[d] - self.lo[d]) > 0.0 && self.bd(d) > 1
    }
    /// block extent along d (1 on inactive axes)
    fn bd(&self, d: usize) -> usize {
        if self.cfg.base_blocks[d] == 0 {
            1
        } else {
            self.cfg.b
        }
    }
    fn dims(&self) -> Vec<usize> {
        (0..3).filter(|&d| self.active(d)).collect()
    }
    /// cells per side at level l
    pub fn n_level(&self, l: usize, d: usize) -> i64 {
        if self.active(d) {
            (self.cfg.base_blocks[d] * self.cfg.b) as i64 * (1 << l)
        } else {
            1
        }
    }
    pub fn dx(&self, l: usize, d: usize) -> f64 {
        (self.hi[d] - self.lo[d]) / self.n_level(l, d) as f64
    }
    fn bcells(&self) -> usize {
        self.bd(0) * self.bd(1) * self.bd(2)
    }
    #[inline]
    fn cidx(&self, c: [usize; 3]) -> usize {
        (c[0] * self.bd(1) + c[1]) * self.bd(2) + c[2]
    }
    /// Cell centre of global cell g at level l.
    pub fn x(&self, l: usize, g: [i64; 3]) -> [f64; 3] {
        std::array::from_fn(|d| if self.active(d) { self.lo[d] + (g[d] as f64 + 0.5) * self.dx(l, d) } else { 0.5 * (self.lo[d] + self.hi[d]) })
    }

    /// Build with `init` (primitive state at a point), refining until the
    /// details of the initial field are resolved to ε.
    pub fn new(sys: S, cfg: AmrConfig, lo: [f64; 3], hi: [f64; 3], init: impl Fn([f64; 3]) -> V<N> + Sync) -> Self {
        let mut a = Amr {
            sys,
            cfg,
            lo,
            hi,
            bc: [[Bc::Periodic; 2]; 3],
            mode: FluxMode::Hybrid,
            shock: ShockScheme::Wcns,
            tol: 1e-3,
            cfl: 0.4,
            positivity: true,
            t: 0.0,
            steps: 0,
            stats: Stats::default(),
            blocks: Vec::new(),
            map: BlockIndex::default(),
            ng: cfg.ghosts.max(5),
            wp: [lagrange6(-3.0, -0.25), lagrange6(-2.0, 0.25)],
            wr: [3.0 / 256.0, -25.0 / 256.0, 150.0 / 256.0, 150.0 / 256.0, -25.0 / 256.0, 3.0 / 256.0],
            stage_dt: 0.0,
            prof: [0.0; 6],
            subcycle: false,
            level_t0: vec![0.0; cfg.levels],
            level_dt: vec![0.0; cfg.levels],
            creg: Default::default(),
            freg: Default::default(),
            regrid_interval: cfg.regrid_every,
            next_regrid: cfg.regrid_every,
            relax: None,
        };
        let nb = cfg.base_blocks.map(|v| v.max(1) as i64);
        a.map = BlockIndex::new(nb, std::array::from_fn(|d| a.active(d)), cfg.levels);
        for i in 0..nb[0] {
            for j in 0..nb[1] {
                for k in 0..nb[2] {
                    a.add_block(0, [i, j, k], &init);
                }
            }
        }
        // refine on the initial data, re-evaluating `init` on new blocks
        for _ in 0..cfg.levels {
            a.restrict_all();
            let flags = a.flags();
            let refine: Vec<usize> = flags.iter().enumerate().filter(|(i, f)| **f > 0 && a.blocks[*i].leaf).map(|(i, _)| i).collect();
            if refine.is_empty() {
                break;
            }
            let keys: Vec<Key> = refine.iter().map(|&i| (a.blocks[i].level as u8, a.blocks[i].idx)).collect();
            for k in keys {
                if let Some(&id) = a.map.get(&k) {
                    a.refine_block(id, Some(&init));
                }
            }
        }
        a
    }

    fn new_block(&self, level: usize, idx: [i64; 3]) -> Block<N> {
        let nc = self.bcells();
        let fsz = |d: usize| -> usize { (0..3).filter(|&e| e != d).map(|e| self.bd(e)).product() };
        Block {
            level,
            idx,
            leaf: true,
            u: vec![[0.0; N]; nc],
            u0: Vec::new(),
            u_save: Vec::new(),
            rhs: vec![[0.0; N]; nc],
            face: std::array::from_fn(|d| [vec![[0.0; N]; fsz(d)], vec![[0.0; N]; fsz(d)]]),
        }
    }

    fn add_block(&mut self, level: usize, idx: [i64; 3], init: &(impl Fn([f64; 3]) -> V<N> + Sync)) {
        let mut b = self.new_block(level, idx);
        let (b0, b1, b2) = (self.bd(0), self.bd(1), self.bd(2));
        for i in 0..b0 {
            for j in 0..b1 {
                for k in 0..b2 {
                    let g = [idx[0] * b0 as i64 + i as i64, idx[1] * b1 as i64 + j as i64, idx[2] * b2 as i64 + k as i64];
                    b.u[(i * b1 + j) * b2 + k] = self.sys.to_cons(&init(self.x(level, g)));
                }
            }
        }
        self.map.insert((level as u8, idx), self.blocks.len());
        self.blocks.push(b);
    }

    #[inline]
    fn block_key(&self, g: [i64; 3]) -> [i64; 3] {
        std::array::from_fn(|d| if self.active(d) { g[d].div_euclid(self.bd(d) as i64) } else { 0 })
    }

    /// Conserved value of global cell g at level l, sampled from the
    /// hierarchy (boundary conditions applied outside the domain).
    pub fn sample(&self, l: usize, g: [i64; 3]) -> V<N> {
        let mut g = g;
        let mut flip = [false; 3];
        // last axis first: a cell outside along several characteristic
        // boundaries then nests them like the axis sweeps of `bc_box`
        // (outermost = last axis)
        for d in (0..3).rev() {
            if !self.active(d) {
                g[d] = 0;
                continue;
            }
            let n = self.n_level(l, d);
            if g[d] < 0 || g[d] >= n {
                let side = (g[d] >= n) as usize;
                match self.bc[d][side] {
                    Bc::Periodic => g[d] = g[d].rem_euclid(n),
                    Bc::Outflow => g[d] = g[d].clamp(0, n - 1),
                    Bc::Reflect => {
                        g[d] = if side == 0 { -1 - g[d] } else { 2 * n - 1 - g[d] };
                        flip[d] = true;
                    }
                    Bc::Dirichlet(v) => return v,
                    Bc::Characteristic(far) => {
                        let mut gi = g;
                        gi[d] = g[d].clamp(0, n - 1);
                        let mut g1 = gi;
                        g1[d] += if side == 1 { -1 } else { 1 };
                        let m = (g[d] - gi[d]).abs() as f64;
                        let (ub, ub1) = (self.sample(l, gi), self.sample(l, g1));
                        let mut v = crate::boundary::characteristic_ghost_layer(&self.sys, &ub, &ub1, &far, d, side, m);
                        // reflections already applied on later axes (the
                        // far state is at rest, so they commute)
                        for e in 0..3 {
                            if flip[e] {
                                v[self.sys.momentum_index(e)] *= -1.0;
                            }
                        }
                        return v;
                    }
                }
            }
        }
        let mut v = match self.map.get(&(l as u8, self.block_key(g))) {
            Some(&id) => {
                let c: [usize; 3] = std::array::from_fn(|d| if self.active(d) { g[d].rem_euclid(self.bd(d) as i64) as usize } else { 0 });
                self.blocks[id].u[self.cidx(c)]
            }
            None => self.predict(l, g),
        };
        for d in 0..3 {
            if flip[d] {
                v[self.sys.momentum_index(d)] *= -1.0;
            }
        }
        v
    }

    /// Fill the cells of a box (level l, global origin g0, size sz) that lie
    /// outside a non-periodic boundary from cells of the same box, one axis
    /// at a time (a cell outside along several axes is filled in the sweep
    /// of its last such axis, from a source already filled in an earlier
    /// sweep). Equivalent to `sample` for those cells but O(1) each: the
    /// per-cell route re-samples the boundary cells through the hierarchy,
    /// which for a characteristic boundary next to a coarser region costs
    /// two 6^D-point predictions per level of recursion. Cells whose source
    /// is not in the box (or not `ok`) are left with ok = false.
    fn bc_box(&self, l: usize, g0: [i64; 3], sz: [usize; 3], val: &mut [V<N>], ok: &mut [bool]) {
        let dims = self.dims();
        let n: [i64; 3] = std::array::from_fn(|d| self.n_level(l, d));
        let at = |c: [i64; 3]| -> Option<usize> {
            for d in 0..3 {
                if c[d] < 0 || c[d] >= sz[d] as i64 {
                    return None;
                }
            }
            Some(((c[0] as usize) * sz[1] + c[1] as usize) * sz[2] + c[2] as usize)
        };
        let outside = |g: [i64; 3], d: usize| -> bool { (g[d] < 0 || g[d] >= n[d]) && self.bc[d][(g[d] >= n[d]) as usize] != Bc::Periodic };
        for (di, &d) in dims.iter().enumerate() {
            for i in 0..sz[0] {
                for j in 0..sz[1] {
                    for k in 0..sz[2] {
                        let c = [i as i64, j as i64, k as i64];
                        let g: [i64; 3] = std::array::from_fn(|e| g0[e] + c[e]);
                        if !outside(g, d) || dims[di + 1..].iter().any(|&e| outside(g, e)) {
                            continue;
                        }
                        let q = at(c).unwrap();
                        let side = (g[d] >= n[d]) as usize;
                        let local = |gg: i64| -> Option<usize> {
                            let mut s = c;
                            s[d] = gg - g0[d];
                            at(s).filter(|&p| ok[p])
                        };
                        let gi = g[d].clamp(0, n[d] - 1);
                        let v = match self.bc[d][side] {
                            Bc::Periodic => unreachable!(),
                            Bc::Dirichlet(v) => Some(v),
                            Bc::Outflow => local(gi).map(|p| val[p]),
                            Bc::Reflect => local(if side == 0 { -1 - g[d] } else { 2 * n[d] - 1 - g[d] }).map(|p| {
                                let mut v = val[p];
                                v[self.sys.momentum_index(d)] *= -1.0;
                                v
                            }),
                            Bc::Characteristic(far) => {
                                let g1 = gi + if side == 1 { -1 } else { 1 };
                                match (local(gi), local(g1)) {
                                    (Some(pb), Some(p1)) => Some(crate::boundary::characteristic_ghost_layer(&self.sys, &val[pb], &val[p1], &far, d, side, (g[d] - gi).abs() as f64)),
                                    _ => None,
                                }
                            }
                        };
                        if let Some(v) = v {
                            val[q] = v;
                            ok[q] = true;
                        }
                    }
                }
            }
        }
    }

    /// `sample` over a box, with the physical boundary cells filled by
    /// `bc_box` (falls back to `sample` where the box lacks the source).
    fn gather_box(&self, l: usize, g0: [i64; 3], sz: [usize; 3]) -> Vec<V<N>> {
        let n: [i64; 3] = std::array::from_fn(|d| self.n_level(l, d));
        let dims = self.dims();
        let len = sz[0] * sz[1] * sz[2];
        let mut val = vec![[0.0; N]; len];
        let mut ok = vec![true; len];
        let mut any = false;
        for i in 0..sz[0] {
            for j in 0..sz[1] {
                for k in 0..sz[2] {
                    let g = [g0[0] + i as i64, g0[1] + j as i64, g0[2] + k as i64];
                    let q = (i * sz[1] + j) * sz[2] + k;
                    if dims.iter().any(|&d| (g[d] < 0 || g[d] >= n[d]) && self.bc[d][(g[d] >= n[d]) as usize] != Bc::Periodic) {
                        ok[q] = false;
                        any = true;
                    } else {
                        val[q] = self.sample(l, g);
                    }
                }
            }
        }
        if any {
            self.bc_box(l, g0, sz, &mut val, &mut ok);
            for i in 0..sz[0] {
                for j in 0..sz[1] {
                    for k in 0..sz[2] {
                        let q = (i * sz[1] + j) * sz[2] + k;
                        if !ok[q] {
                            val[q] = self.sample(l, [g0[0] + i as i64, g0[1] + j as i64, g0[2] + k as i64]);
                        }
                    }
                }
            }
        }
        val
    }

    /// 6th-order prediction of the level-l box [fmin, fmax] from level
    /// l − 1: one gathered coarse patch, one 6-point pass per active axis.
    /// Returns the values (row-major) and the box extents.
    fn predict_box(&self, l: usize, fmin: [i64; 3], fmax: [i64; 3]) -> (Vec<V<N>>, [usize; 3]) {
        let dims = self.dims();
        let act: [bool; 3] = std::array::from_fn(|d| dims.contains(&d));
        let cmin: [i64; 3] = std::array::from_fn(|d| if act[d] { fmin[d].div_euclid(2) - 3 } else { 0 });
        let cmax: [i64; 3] = std::array::from_fn(|d| if act[d] { fmax[d].div_euclid(2) + 3 } else { 0 });
        let mut sz: [usize; 3] = std::array::from_fn(|d| (cmax[d] - cmin[d] + 1) as usize);
        let mut arr = self.gather_box(l - 1, cmin, sz);
        // passes along each active dimension: coarse index → fine index
        for &d in &dims {
            let nf = (fmax[d] - fmin[d] + 1) as usize;
            let mut nsz = sz;
            nsz[d] = nf;
            let mut na = vec![[0.0; N]; nsz[0] * nsz[1] * nsz[2]];
            for i in 0..nsz[0] {
                for j in 0..nsz[1] {
                    for k in 0..nsz[2] {
                        let c = [i, j, k];
                        let f = fmin[d] + c[d] as i64;
                        let par = f.rem_euclid(2) as usize;
                        let base = f.div_euclid(2) + if par == 1 { -2 } else { -3 } - cmin[d];
                        let mut acc = [0.0; N];
                        for kk in 0..6 {
                            let mut src = c;
                            src[d] = (base + kk as i64) as usize;
                            let v = &arr[(src[0] * sz[1] + src[1]) * sz[2] + src[2]];
                            let wgt = self.wp[par][kk];
                            for q in 0..N {
                                acc[q] += wgt * v[q];
                            }
                        }
                        na[(i * nsz[1] + j) * nsz[2] + k] = acc;
                    }
                }
            }
            arr = na;
            sz = nsz;
        }
        (arr, sz)
    }

    /// 6th-order prediction of cell g (level l) from level l − 1.
    fn predict(&self, l: usize, g: [i64; 3]) -> V<N> {
        let dims = self.dims();
        let mut acc = [0.0; N];
        let nd = dims.len();
        let total = 6usize.pow(nd as u32);
        for s in 0..total {
            let mut w = 1.0;
            let mut c = [0i64; 3];
            let mut rem = s;
            for &d in &dims {
                let k = rem % 6;
                rem /= 6;
                let odd = g[d].rem_euclid(2) as usize;
                let base = g[d].div_euclid(2) + if odd == 1 { -2 } else { -3 };
                c[d] = base + k as i64;
                w *= self.wp[odd][k];
            }
            let v = self.sample(l - 1, c);
            for q in 0..N {
                acc[q] += w * v[q];
            }
        }
        acc
    }

    /// Separable 6th-order restriction of ancestor block (l, idx) from level
    /// l + 1: one fine patch is gathered (with a validity mask: present at
    /// l + 1, or supplied by a non-periodic boundary condition), then one
    /// pass per direction maps fine to coarse indices with the centred
    /// 6-point midpoint stencil, shifted (o = ±1, ±2) along that line where
    /// the fine data ends; the child average is the last resort.
    /// Computes the local cell box [lo, hi) only; returns the box's values in
    /// row-major order.
    fn restrict_box(&self, l: usize, idx: [i64; 3], lo: [usize; 3], hi: [usize; 3]) -> Vec<V<N>> {
        let dims = self.dims();
        let bd0 = [self.bd(0), self.bd(1), self.bd(2)];
        let bd: [usize; 3] = std::array::from_fn(|d| hi[d] - lo[d]);
        let act: [bool; 3] = std::array::from_fn(|d| dims.contains(&d));
        let c0: [i64; 3] = std::array::from_fn(|d| idx[d] * bd0[d] as i64 + lo[d] as i64);
        let f0: [i64; 3] = std::array::from_fn(|d| if act[d] { 2 * c0[d] - 4 } else { 0 });
        let mut sz: [usize; 3] = std::array::from_fn(|d| if act[d] { 2 * bd[d] + 8 } else { 1 });
        let nf: [i64; 3] = std::array::from_fn(|d| self.n_level(l + 1, d));
        let mut val: Vec<V<N>> = Vec::with_capacity(sz[0] * sz[1] * sz[2]);
        let mut ok: Vec<bool> = Vec::with_capacity(sz[0] * sz[1] * sz[2]);
        let mut any_out = false;
        for i in 0..sz[0] {
            for j in 0..sz[1] {
                for k in 0..sz[2] {
                    let g = [f0[0] + i as i64, f0[1] + j as i64, f0[2] + k as i64];
                    let mut gw = g;
                    let mut outside = false;
                    for d in 0..3 {
                        if act[d] && (g[d] < 0 || g[d] >= nf[d]) {
                            if self.bc[d][(g[d] >= nf[d]) as usize] == Bc::Periodic {
                                gw[d] = g[d].rem_euclid(nf[d]);
                            } else {
                                outside = true;
                            }
                        }
                    }
                    if outside {
                        val.push([0.0; N]);
                        ok.push(false);
                        any_out = true;
                    } else if let Some(&bid) = self.map.get(&((l + 1) as u8, self.block_key(gw))) {
                        let c: [usize; 3] = std::array::from_fn(|d| if act[d] { gw[d].rem_euclid(bd0[d] as i64) as usize } else { 0 });
                        val.push(self.blocks[bid].u[self.cidx(c)]);
                        ok.push(true);
                    } else {
                        val.push([0.0; N]);
                        ok.push(false);
                    }
                }
            }
        }
        if any_out {
            // boundary cells whose sources are not refined stay invalid (the
            // stencil shifts inwards) rather than being predicted from the
            // coarser level through the boundary condition
            self.bc_box(l + 1, f0, sz, &mut val, &mut ok);
        }
        let shifts: [i64; 5] = [0, -1, 1, -2, 2];
        let wts: [[f64; 6]; 5] = std::array::from_fn(|m| if m == 0 { self.wr } else { lagrange6(-2.5 + shifts[m] as f64, 0.0) });
        let mut off = f0;
        for &d in &dims {
            let mut nsz = sz;
            nsz[d] = bd[d];
            let mut nv = vec![[0.0; N]; nsz[0] * nsz[1] * nsz[2]];
            let mut nok = vec![false; nsz[0] * nsz[1] * nsz[2]];
            for i in 0..nsz[0] {
                for j in 0..nsz[1] {
                    for k in 0..nsz[2] {
                        let c = [i, j, k];
                        let at_src = |pos: usize| -> usize {
                            let mut s2 = c;
                            s2[d] = pos;
                            (s2[0] * sz[1] + s2[1]) * sz[2] + s2[2]
                        };
                        // local fine index of the centred stencil's first node
                        let base = 2 * c[d] as i64 + 2; // (2(c0+c) − 2) − (2c0 − 4)
                        let mut done = false;
                        for (m, &o) in shifts.iter().enumerate() {
                            let st = base + o;
                            if st < 0 || st + 5 >= sz[d] as i64 {
                                continue;
                            }
                            if (0..6).all(|kk| ok[at_src((st + kk) as usize)]) {
                                let mut acc = [0.0; N];
                                for kk in 0..6 {
                                    let v = &val[at_src((st + kk) as usize)];
                                    for q in 0..N {
                                        acc[q] += wts[m][kk as usize] * v[q];
                                    }
                                }
                                nv[(i * nsz[1] + j) * nsz[2] + k] = acc;
                                nok[(i * nsz[1] + j) * nsz[2] + k] = true;
                                done = true;
                                break;
                            }
                        }
                        if !done {
                            let (pa, pb) = (at_src((base + 2) as usize), at_src((base + 3) as usize));
                            if ok[pa] && ok[pb] {
                                let mut acc = [0.0; N];
                                for q in 0..N {
                                    acc[q] = 0.5 * (val[pa][q] + val[pb][q]);
                                }
                                nv[(i * nsz[1] + j) * nsz[2] + k] = acc;
                                nok[(i * nsz[1] + j) * nsz[2] + k] = true;
                            }
                        }
                    }
                }
            }
            val = nv;
            ok = nok;
            sz = nsz;
            off[d] = c0[d];
        }
        let _ = off;
        val
    }

    /// Refresh the ancestor blocks that border coarser leaves ("edge"
    /// ancestors), finest first. Only those are ever sampled — by ghost
    /// cells next to the coarse region and by predictions into it — and by
    /// 2:1 balance their restriction stencils near that edge read leaves,
    /// so the cost scales with the surface of the refined region.
    fn restrict_all(&mut self) {
        for l in (0..self.cfg.levels.saturating_sub(1)).rev() {
            let ids: Vec<usize> = (0..self.blocks.len())
                .filter(|&i| {
                    let b = &self.blocks[i];
                    !b.leaf && b.level == l && self.neighbours(l, b.idx).iter().any(|nk| self.map.get(&(l as u8, *nk)).map_or(true, |&j| self.blocks[j].leaf))
                })
                .collect();
            // only the cells within the ghost depth of a non-ancestor
            // neighbour are ever read (ghost copies into same-level leaves,
            // predictions into finer blocks): restrict those boxes only
            let ng = self.ng;
            let dims = self.dims();
            let bd = [self.bd(0), self.bd(1), self.bd(2)];
            let boxes: Vec<(usize, Vec<([usize; 3], [usize; 3])>)> = ids
                .iter()
                .map(|&id| {
                    let b = &self.blocks[id];
                    let nbk: [i64; 3] = std::array::from_fn(|d| if dims.contains(&d) { self.cfg.base_blocks[d] as i64 * (1 << l) } else { 1 });
                    let mut bx = Vec::new();
                    let offs = |d: usize| -> Vec<i64> { if dims.contains(&d) { vec![-1, 0, 1] } else { vec![0] } };
                    for &o0 in &offs(0) {
                        for &o1 in &offs(1) {
                            for &o2 in &offs(2) {
                                let o = [o0, o1, o2];
                                if o == [0, 0, 0] {
                                    continue;
                                }
                                let mut nk = [0i64; 3];
                                let mut skip = false;
                                for d in 0..3 {
                                    nk[d] = b.idx[d] + o[d];
                                    if nk[d] < 0 || nk[d] >= nbk[d] {
                                        if self.bc[d][(nk[d] >= nbk[d]) as usize] == Bc::Periodic {
                                            nk[d] = nk[d].rem_euclid(nbk[d]);
                                        } else {
                                            skip = true;
                                        }
                                    }
                                }
                                if skip || self.map.get(&(l as u8, nk)).is_some_and(|&j| !self.blocks[j].leaf) {
                                    continue;
                                }
                                let w = (ng + 3).min(bd[0].max(bd[1]).max(bd[2]));
                                let lo: [usize; 3] = std::array::from_fn(|d| if o[d] == 1 { bd[d].saturating_sub(w) } else { 0 });
                                let hi: [usize; 3] = std::array::from_fn(|d| if o[d] == -1 { w.min(bd[d]) } else { bd[d] });
                                bx.push((lo, hi));
                            }
                        }
                    }
                    // many overlapping slabs cost more than one whole-block pass
                    let vol: usize = bx.iter().map(|(lo, hi)| (0..3).map(|d| hi[d] - lo[d]).product::<usize>()).sum();
                    if 2 * vol > bd[0] * bd[1] * bd[2] {
                        bx = vec![([0; 3], bd)];
                    }
                    (id, bx)
                })
                .collect();
            let new: Vec<(usize, Vec<(([usize; 3], [usize; 3]), Vec<V<N>>)>)> = boxes
                .par_iter()
                .map(|(id, bx)| {
                    let idx = self.blocks[*id].idx;
                    (*id, bx.iter().map(|&(lo, hi)| ((lo, hi), self.restrict_box(l, idx, lo, hi))).collect())
                })
                .collect();
            for (id, parts) in new {
                for ((lo, hi), vals) in parts {
                    let mut it = vals.into_iter();
                    for i in lo[0]..hi[0] {
                        for j in lo[1]..hi[1] {
                            for k in lo[2]..hi[2] {
                                let c = self.cidx([i, j, k]);
                                self.blocks[id].u[c] = it.next().unwrap();
                            }
                        }
                    }
                }
            }
        }
    }

    /// Largest relative wavelet detail of a block's cells: in each active
    /// direction, the difference between u_i and its 6th-order midpoint
    /// prediction from the neighbours at spacing 2h (i ± 1, ± 3, ± 5), i.e.
    /// the error of representing the data one level coarser. Measured on
    /// the sensor values (relative).
    fn block_detail(&self, id: usize) -> f64 {
        let dims = self.dims();
        let ng = self.ng;
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let act: [bool; 3] = std::array::from_fn(|d| dims.contains(&d));
        let ext: [usize; 3] = std::array::from_fn(|d| if act[d] { bd[d] + 2 * ng } else { 1 });
        let gof: [usize; 3] = std::array::from_fn(|d| if act[d] { ng } else { 0 });
        let scr = self.fill_scratch(id);
        // sensor values on the scratch box (needs ±5 cells: ng ≥ 5)
        let mut ns = 0;
        let sens: Vec<[f64; 8]> = scr
            .iter()
            .map(|u| {
                let mut v = [0.0; 8];
                ns = self.sys.sensor(&self.sys.to_prim(u), &mut v);
                v
            })
            .collect();
        let at = |c: [usize; 3]| (c[0] * ext[1] + c[1]) * ext[2] + c[2];
        let mut mx: f64 = 0.0;
        for i in 0..bd[0] {
            for j in 0..bd[1] {
                for k in 0..bd[2] {
                    let c = [i + gof[0], j + gof[1], k + gof[2]];
                    let s0 = &sens[at(c)];
                    for &d in &dims {
                        let mut pred = [0.0; 8];
                        for (m, off) in [-5i64, -3, -1, 1, 3, 5].iter().enumerate() {
                            let mut cc = c;
                            cc[d] = (cc[d] as i64 + off) as usize;
                            let sv = &sens[at(cc)];
                            for q in 0..ns {
                                pred[q] += self.wr[m] * sv[q];
                            }
                        }
                        for q in 0..ns {
                            let dd = (s0[q] - pred[q]).abs();
                            mx = if dd.is_finite() { mx.max(dd) } else { f64::INFINITY };
                        }
                    }
                }
            }
        }
        mx
    }

    /// Refinement flags: +1 refine, −1 may coarsen, 0 keep.
    fn flags(&self) -> Vec<i8> {
        let det: Vec<f64> = (0..self.blocks.len()).into_par_iter().map(|i| if self.blocks[i].leaf { self.block_detail(i) } else { 0.0 }).collect();
        self.flags_from_details(&det)
    }

    /// Flags from per-block details (block id order; ancestors ignored),
    /// e.g. computed on an accelerator.
    pub fn flags_from_details(&self, det: &[f64]) -> Vec<i8> {
        (0..self.blocks.len())
            .map(|i| {
                let b = &self.blocks[i];
                if !b.leaf {
                    0
                } else if det[i] > self.cfg.eps && b.level + 1 < self.cfg.levels {
                    1
                } else if det[i] < self.cfg.eps / 128.0 {
                    -1
                } else {
                    0
                }
            })
            .collect()
    }

    fn children(&self, idx: [i64; 3]) -> Vec<[i64; 3]> {
        let dims = self.dims();
        let mut out = vec![[0i64; 3]];
        for d in 0..3 {
            let mut next = Vec::new();
            for c in &out {
                if dims.contains(&d) {
                    for o in 0..2 {
                        let mut cc = *c;
                        cc[d] = 2 * idx[d] + o;
                        next.push(cc);
                    }
                } else {
                    let mut cc = *c;
                    cc[d] = 0;
                    next.push(cc);
                }
            }
            out = next;
        }
        out
    }

    /// Neighbouring block keys (26 in 3D) at the same level, periodic wrap.
    fn neighbours(&self, l: usize, idx: [i64; 3]) -> Vec<[i64; 3]> {
        let dims = self.dims();
        let nb: [i64; 3] = std::array::from_fn(|d| if self.active(d) { self.cfg.base_blocks[d] as i64 * (1 << l) } else { 1 });
        let mut out = Vec::new();
        let r = |d: usize| if dims.contains(&d) { -1..=1 } else { 0..=0 };
        for a in r(0) {
            for b in r(1) {
                for c in r(2) {
                    if a == 0 && b == 0 && c == 0 {
                        continue;
                    }
                    let mut k = [idx[0] + a, idx[1] + b, idx[2] + c];
                    let mut inside = true;
                    for d in 0..3 {
                        if k[d] < 0 || k[d] >= nb[d] {
                            if self.bc[d][(k[d] >= nb[d]) as usize] == Bc::Periodic {
                                k[d] = k[d].rem_euclid(nb[d]);
                            } else {
                                inside = false;
                            }
                        }
                    }
                    if inside {
                        out.push(k);
                    }
                }
            }
        }
        out
    }

    /// Refine leaf `id` into 2^D children (data predicted, or re-evaluated
    /// from `init` when given), enforcing 2:1 balance first. The parent
    /// stays as an ancestor.
    fn refine_block(&mut self, id: usize, init: Option<&(dyn Fn([f64; 3]) -> V<N> + Sync)>) {
        let (l, idx) = (self.blocks[id].level, self.blocks[id].idx);
        if !self.blocks[id].leaf {
            return;
        }
        // balance: every same-level neighbour position must exist (a block
        // at level l, leaf or ancestor); otherwise refine the coarser leaf
        // covering it first
        if l > 0 {
            for nk in self.neighbours(l, idx) {
                if !self.map.contains_key(&(l as u8, nk)) {
                    // the covering leaf is at a coarser level m < l
                    for m in (0..l).rev() {
                        let pk: [i64; 3] = std::array::from_fn(|d| if self.active(d) { nk[d] >> (l - m) } else { 0 });
                        if let Some(&pid) = self.map.get(&(m as u8, pk)) {
                            self.refine_block(pid, init);
                            break;
                        }
                    }
                }
            }
        }
        let Some(&id) = self.map.get(&(l as u8, idx)) else { return };
        if !self.blocks[id].leaf {
            return;
        }
        let kids = self.children(idx);
        let (b0, b1, b2) = (self.bd(0), self.bd(1), self.bd(2));
        let bd = [b0, b1, b2];
        let newb: Vec<Block<N>> = kids
            .par_iter()
            .map(|&c| {
                let mut b = self.new_block(l + 1, c);
                match init {
                    Some(f) => {
                        for i in 0..b0 {
                            for j in 0..b1 {
                                for k in 0..b2 {
                                    let g = [c[0] * b0 as i64 + i as i64, c[1] * b1 as i64 + j as i64, c[2] * b2 as i64 + k as i64];
                                    b.u[(i * b1 + j) * b2 + k] = self.sys.to_cons(&f(self.x(l + 1, g)));
                                }
                            }
                        }
                    }
                    None => {
                        // separable box prediction (same operator as `predict`)
                        let fmin: [i64; 3] = std::array::from_fn(|d| c[d] * bd[d] as i64);
                        let fmax: [i64; 3] = std::array::from_fn(|d| c[d] * bd[d] as i64 + bd[d] as i64 - 1);
                        b.u = self.predict_box(l + 1, fmin, fmax).0;
                    }
                }
                b
            })
            .collect();
        self.blocks[id].leaf = false;
        for b in newb {
            self.map.insert((b.level as u8, b.idx), self.blocks.len());
            self.blocks.push(b);
        }
    }

    /// Merge the 2^D leaf children of ancestor (l, idx) back into it.
    fn coarsen(&mut self, l: usize, idx: [i64; 3]) {
        let kids = self.children(idx);
        let pid = self.map[&(l as u8, idx)];
        let (b0, b1, b2) = (self.bd(0), self.bd(1), self.bd(2));
        let u = self.restrict_box(l, idx, [0; 3], [b0, b1, b2]);
        self.blocks[pid].u = u;
        self.blocks[pid].leaf = true;
        let mut ids: Vec<usize> = kids.iter().map(|c| self.map[&((l + 1) as u8, *c)]).collect();
        ids.sort_unstable_by(|a, b| b.cmp(a));
        for i in ids {
            let b = &self.blocks[i];
            self.map.remove(&(b.level as u8, b.idx));
            self.blocks.swap_remove(i);
            if i < self.blocks.len() {
                let m = &self.blocks[i];
                self.map.insert((m.level as u8, m.idx), i);
            }
        }
    }

    fn leaves(&self) -> Vec<usize> {
        (0..self.blocks.len()).filter(|&i| self.blocks[i].leaf).collect()
    }

    /// Adapt the hierarchy to the current solution.
    pub fn regrid(&mut self) {
        self.restrict_all();
        let flags = self.flags();
        let protected = self.regrid_refine(&flags);
        self.restrict_all();
        let flags = self.flags();
        self.regrid_coarsen(&flags, &protected);
    }

    /// First half of `regrid`: refine flagged leaves and their same-level
    /// leaf neighbours (buffer). Returns the blocks protected from
    /// coarsening in the second half. Ancestors must be restricted.
    pub fn regrid_refine(&mut self, flags: &[i8]) -> std::collections::HashSet<(u8, [i64; 3])> {
        // refine flagged leaves and their same-level neighbours (buffer)
        let mut want: Vec<Key> = Vec::new();
        for (i, f) in flags.iter().enumerate() {
            if *f > 0 {
                let b = &self.blocks[i];
                want.push((b.level as u8, b.idx));
                for nk in self.neighbours(b.level, b.idx) {
                    if self.map.get(&(b.level as u8, nk)).is_some_and(|&j| self.blocks[j].leaf) {
                        want.push((b.level as u8, nk));
                    }
                }
            }
        }
        want.sort();
        want.dedup();
        let protected: std::collections::HashSet<Key> = want.iter().cloned().collect();
        for k in want {
            if let Some(&id) = self.map.get(&k) {
                if (k.0 as usize) + 1 < self.cfg.levels {
                    self.refine_block(id, None);
                }
            }
        }
        protected
    }

    /// Second half of `regrid`: coarsen ancestors whose children are all
    /// leaves with tiny details (`flags` of the refined hierarchy).
    pub fn regrid_coarsen(&mut self, flags: &[i8], protected: &std::collections::HashSet<(u8, [i64; 3])>) {
        let nkids = 1usize << self.dims().len();
        let mut cand: Vec<Key> = Vec::new();
        for b in self.blocks.iter().filter(|b| !b.leaf) {
            let kids = self.children(b.idx);
            let ok = kids.len() == nkids
                && kids.iter().all(|c| match self.map.get(&((b.level + 1) as u8, *c)) {
                    Some(&j) => self.blocks[j].leaf && flags[j] < 0 && !protected.contains(&((b.level + 1) as u8, *c)),
                    None => false,
                });
            if ok {
                cand.push((b.level as u8, b.idx));
            }
        }
        cand.sort_by(|a, b| b.0.cmp(&a.0));
        for (l, idx) in cand {
            let l = l as usize;
            // keep 2:1 balance: no neighbour of the parent may hold blocks at
            // level l + 2 (i.e. neighbours' children must all be leaves)
            let ok = self.neighbours(l, idx).iter().all(|nk| {
                self.children(*nk).iter().all(|c| self.map.get(&((l + 1) as u8, *c)).map_or(true, |&j| self.blocks[j].leaf))
            });
            if ok && self.map.contains_key(&(l as u8, idx)) {
                self.coarsen(l, idx);
            }
        }
    }

    /// Conserved values of leaf `id` on its ghost-extended box, filled block
    /// by block: a neighbour block at the same level (leaf or ancestor) is a
    /// slab copy; a coarser region is predicted from one gathered coarse
    /// patch by separable 6th-order interpolation (3 passes of 6 points
    /// instead of 6^D per cell); boundary slabs fall back to `sample`.
    fn fill_scratch(&self, id: usize) -> Vec<V<N>> {
        let b = &self.blocks[id];
        let l = b.level;
        let ng = self.ng;
        let dims = self.dims();
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let act: [bool; 3] = std::array::from_fn(|d| dims.contains(&d));
        let ext: [usize; 3] = std::array::from_fn(|d| if act[d] { bd[d] + 2 * ng } else { 1 });
        let gof: [usize; 3] = std::array::from_fn(|d| if act[d] { ng } else { 0 });
        let at = |c: [usize; 3]| (c[0] * ext[1] + c[1]) * ext[2] + c[2];
        let mut out = vec![[0.0; N]; ext[0] * ext[1] * ext[2]];
        let mut filled = vec![true; out.len()];
        for i in 0..bd[0] {
            for j in 0..bd[1] {
                for k in 0..bd[2] {
                    out[at([i + gof[0], j + gof[1], k + gof[2]])] = b.u[self.cidx([i, j, k])];
                }
            }
        }
        let nbk: [i64; 3] = std::array::from_fn(|d| if act[d] { self.cfg.base_blocks[d] as i64 * (1 << l) } else { 1 });
        let rng = |o: i64, d: usize| -> std::ops::Range<usize> {
            if !act[d] {
                return 0..1;
            }
            match o {
                -1 => 0..ng,
                0 => ng..ng + bd[d],
                _ => ng + bd[d]..2 * ng + bd[d],
            }
        };
        let offs = |d: usize| -> Vec<i64> { if act[d] { vec![-1, 0, 1] } else { vec![0] } };
        for &o0 in &offs(0) {
            for &o1 in &offs(1) {
                for &o2 in &offs(2) {
                    let o = [o0, o1, o2];
                    if o == [0, 0, 0] {
                        continue;
                    }
                    let r: [std::ops::Range<usize>; 3] = [rng(o0, 0), rng(o1, 1), rng(o2, 2)];
                    let gl = |c: [usize; 3]| -> [i64; 3] { std::array::from_fn(|d| b.idx[d] * bd[d] as i64 + c[d] as i64 - gof[d] as i64) };
                    // neighbour key, periodic wrap; None if outside a non-periodic boundary
                    let mut nk = [0i64; 3];
                    let mut inside = true;
                    for d in 0..3 {
                        nk[d] = b.idx[d] + o[d];
                        if nk[d] < 0 || nk[d] >= nbk[d] {
                            if self.bc[d][(nk[d] >= nbk[d]) as usize] == Bc::Periodic {
                                nk[d] = nk[d].rem_euclid(nbk[d]);
                            } else {
                                inside = false;
                            }
                        }
                    }
                    if !inside {
                        // filled from the box below (physical boundary)
                        for i in r[0].clone() {
                            for j in r[1].clone() {
                                for k in r[2].clone() {
                                    filled[at([i, j, k])] = false;
                                }
                            }
                        }
                        continue;
                    }
                    if let Some(&nid) = self.map.get(&(l as u8, nk)) {
                        let nb = &self.blocks[nid];
                        for i in r[0].clone() {
                            for j in r[1].clone() {
                                for k in r[2].clone() {
                                    let c = [i, j, k];
                                    let lc: [usize; 3] = std::array::from_fn(|d| if act[d] { (c[d] as i64 - gof[d] as i64 - o[d] * bd[d] as i64) as usize } else { 0 });
                                    out[at(c)] = nb.u[self.cidx(lc)];
                                }
                            }
                        }
                        continue;
                    }
                    // coarser region: gather coarse patch, separable prediction
                    let fmin: [i64; 3] = std::array::from_fn(|d| gl([r[0].start, r[1].start, r[2].start])[d]);
                    let fmax: [i64; 3] = std::array::from_fn(|d| gl([r[0].end - 1, r[1].end - 1, r[2].end - 1])[d]);
                    let (arr, sz) = self.predict_box(l, fmin, fmax);
                    for i in r[0].clone() {
                        for j in r[1].clone() {
                            for k in r[2].clone() {
                                let c = [i, j, k];
                                let g = gl(c);
                                let li: [usize; 3] = std::array::from_fn(|d| if act[d] { (g[d] - fmin[d]) as usize } else { 0 });
                                out[at(c)] = arr[(li[0] * sz[1] + li[1]) * sz[2] + li[2]];
                            }
                        }
                    }
                }
            }
        }
        if filled.iter().any(|&f| !f) {
            let g0: [i64; 3] = std::array::from_fn(|d| b.idx[d] * bd[d] as i64 - gof[d] as i64);
            self.bc_box(l, g0, ext, &mut out, &mut filled);
            for i in 0..ext[0] {
                for j in 0..ext[1] {
                    for k in 0..ext[2] {
                        if !filled[at([i, j, k])] {
                            out[at([i, j, k])] = self.sample(l, [g0[0] + i as i64, g0[1] + j as i64, g0[2] + k as i64]);
                        }
                    }
                }
            }
        }
        out
    }

    /// Right-hand side and boundary face fluxes of one block.
    fn block_rhs(&self, id: usize, ghost_ns: &std::sync::atomic::AtomicU64, flux_ns: &std::sync::atomic::AtomicU64) -> (Vec<V<N>>, [[Vec<V<N>>; 2]; 3], Stats) {
        let ng = self.ng;
        let dims = self.dims();
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let ext: [usize; 3] = std::array::from_fn(|d| if dims.contains(&d) { bd[d] + 2 * ng } else { 1 });
        let gof: [usize; 3] = std::array::from_fn(|d| if dims.contains(&d) { ng } else { 0 });
        let ndim = dims.len() as f64;
        let b = &self.blocks[id];
        let l = b.level;
        let tg = std::time::Instant::now();
        // scratch field with ghosts (block-level exchange), then primitives
        let w: Vec<V<N>> = self.fill_scratch(id).iter().map(|u| self.sys.to_prim(u)).collect();
        ghost_ns.fetch_add(tg.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        let tf = std::time::Instant::now();
        let mut rhs = vec![[0.0; N]; self.bcells()];
        let mut face: [[Vec<V<N>>; 2]; 3] = std::array::from_fn(|_| [Vec::new(), Vec::new()]);
        let mut st = Stats::default();
        for &d in &dims {
            let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
            let (e0, e1) = (others[0], others[1]);
            let dx = self.dx(l, d);
            let opts = SchemeOpts { mode: self.mode, shock: self.shock, tol: self.tol, lambda: if self.positivity { ndim * self.stage_dt / dx } else { 0.0 } };
            let nc = self.sys.nc_vars();
            let ncg = self.sys.nc_general();
            let mut fl = vec![[0.0; N]; bd[d] + 1];
            let mut uf = vec![0.0; bd[d] + 1];
            let mut wf = vec![[0.0; N]; bd[d] + 1];
            let mut line = vec![[0.0; N]; ext[d]];
            let mut f_lo = Vec::with_capacity(bd[e0] * bd[e1]);
            let mut f_hi = Vec::with_capacity(bd[e0] * bd[e1]);
            for a in 0..bd[e0] {
                for bb in 0..bd[e1] {
                    let at = |cc: usize| {
                let mut c = [0usize; 3];
                c[d] = cc;
                c[e0] = a + gof[e0];
                c[e1] = bb + gof[e1];
                (c[0] * ext[1] + c[1]) * ext[2] + c[2]
                    };
                    for cc in 0..ext[d] {
                line[cc] = w[at(cc)];
                    }
                    line_fluxes(&self.sys, &line, d, ng, opts, &mut fl, &mut uf, &mut wf, &mut st);
                    for cc in 0..bd[d] {
                let mut c = [0usize; 3];
                c[d] = cc;
                c[e0] = a;
                c[e1] = bb;
                let o = &mut rhs[self.cidx(c)];
                for q in 0..N {
                    o[q] -= (fl[cc + 1][q] - fl[cc][q]) / dx;
                }
                let wc = &line[cc + ng];
                for qq in nc.clone() {
                    o[qq] += self.sys.nc_coef_k(wc, qq) * (uf[cc + 1] - uf[cc]) / dx;
                }
                if ncg {
                    let mut dw = [0.0; N];
                    for q in 0..N {
                        dw[q] = (wf[cc + 1][q] - wf[cc][q]) / dx;
                    }
                    let add = self.sys.nc_rhs(wc, &dw, d);
                    for q in 0..N {
                        o[q] += add[q];
                    }
                }
                    }
                    f_lo.push(fl[0]);
                    f_hi.push(fl[bd[d]]);
                }
            }
            face[d] = [f_lo, f_hi];
        }
        flux_ns.fetch_add(tf.elapsed().as_nanos() as u64, std::sync::atomic::Ordering::Relaxed);
        (rhs, face, st)
    }

    /// Right-hand side of every block, with refluxing at level interfaces.
    fn compute_rhs(&mut self) -> Stats {
        let t0 = std::time::Instant::now();
        self.restrict_all();
        self.prof[0] += t0.elapsed().as_secs_f64();
        let ghost_ns = std::sync::atomic::AtomicU64::new(0);
        let flux_ns = std::sync::atomic::AtomicU64::new(0);
        let leaves = self.leaves();
        let results: Vec<(Vec<V<N>>, [[Vec<V<N>>; 2]; 3], Stats)> = leaves.par_iter().map(|&id| self.block_rhs(id, &ghost_ns, &flux_ns)).collect();
        let nthreads = rayon::current_num_threads() as f64;
        self.prof[1] += ghost_ns.into_inner() as f64 * 1e-9 / nthreads;
        self.prof[2] += flux_ns.into_inner() as f64 * 1e-9 / nthreads;
        let mut total = Stats::default();
        for (&id, (rhs, face, st)) in leaves.iter().zip(results) {
            self.blocks[id].rhs = rhs;
            self.blocks[id].face = face;
            total.add(&st);
        }
        let tr = std::time::Instant::now();
        self.reflux();
        self.prof[3] += tr.elapsed().as_secs_f64();
        total
    }

    /// Replace coarse boundary fluxes by the mean of the finer neighbour's
    /// fluxes over the same face (exact conservation across levels).
    fn reflux(&mut self) {
        let dims = self.dims();
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let mut corr: Vec<(usize, usize, V<N>)> = Vec::new(); // (block, cell, Δrhs)
        for (id, b) in self.blocks.iter().enumerate() {
            let l = b.level;
            if !b.leaf || l + 1 >= self.cfg.levels {
                continue;
            }
            for &d in &dims {
                let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
                let (e0, e1) = (others[0], others[1]);
                let dx = self.dx(l, d);
                let nbk = self.cfg.base_blocks[d] as i64 * (1 << l);
                for side in 0..2 {
                    let mut nk = b.idx;
                    nk[d] += if side == 1 { 1 } else { -1 };
                    if nk[d] < 0 || nk[d] >= nbk {
                        if self.bc[d][side] != Bc::Periodic {
                            continue;
                        }
                        nk[d] = nk[d].rem_euclid(nbk);
                    }
                    // the neighbour is refined iff it is an ancestor
                    if !self.map.get(&(l as u8, nk)).is_some_and(|&j| !self.blocks[j].leaf) {
                        continue;
                    }
                    for a in 0..bd[e0] {
                        for bb in 0..bd[e1] {
                            // fine faces covering coarse face (a, bb)
                            let mut sum = [0.0; N];
                            let mut cnt = 0.0;
                            let fa_r = if dims.contains(&e0) { 0..2 } else { 0..1 };
                            for fa in fa_r {
                                let fb_r = if dims.contains(&e1) { 0..2 } else { 0..1 };
                                for fb in fb_r {
                                    let ga = 2 * (b.idx[e0] * bd[e0] as i64 + a as i64) + fa;
                                    let gb = 2 * (b.idx[e1] * bd[e1] as i64 + bb as i64) + fb;
                                    let mut ck = nk;
                                    ck[d] = 2 * nk[d] + if side == 1 { 0 } else { 1 };
                                    ck[e0] = if dims.contains(&e0) { ga.div_euclid(bd[e0] as i64) } else { 0 };
                                    ck[e1] = if dims.contains(&e1) { gb.div_euclid(bd[e1] as i64) } else { 0 };
                                    let Some(&cid) = self.map.get(&((l + 1) as u8, ck)) else { continue };
                                    if !self.blocks[cid].leaf {
                                        continue;
                                    }
                                    let la = if dims.contains(&e0) { ga.rem_euclid(bd[e0] as i64) as usize } else { 0 };
                                    let lb = if dims.contains(&e1) { gb.rem_euclid(bd[e1] as i64) as usize } else { 0 };
                                    let fv = self.blocks[cid].face[d][1 - side][la * bd[e1] + lb];
                                    for q in 0..N {
                                        sum[q] += fv[q];
                                    }
                                    cnt += 1.0;
                                }
                            }
                            if cnt == 0.0 {
                                continue;
                            }
                            let fc = b.face[d][side][a * bd[e1] + bb];
                            let mut c = [0usize; 3];
                            c[d] = if side == 1 { bd[d] - 1 } else { 0 };
                            c[e0] = a;
                            c[e1] = bb;
                            let mut dr = [0.0; N];
                            for q in 0..N {
                                let df = sum[q] / cnt - fc[q];
                                dr[q] = if side == 1 { -df / dx } else { df / dx };
                            }
                            corr.push((id, self.cidx(c), dr));
                        }
                    }
                }
            }
        }
        for (id, c, dr) in corr {
            for q in 0..N {
                self.blocks[id].rhs[c][q] += dr[q];
            }
        }
    }

    /// Stable global time step (finest cells dominate).
    pub fn dt_cfl(&self) -> f64 {
        let dims = self.dims();
        let rate = self
            .blocks
            .par_iter()
            .filter(|b| b.leaf)
            .map(|b| {
                let mut r: f64 = 0.0;
                for u in &b.u {
                    let w = self.sys.to_prim(u);
                    let s: f64 = dims.iter().map(|&d| self.sys.max_speed(&w, d) / self.dx(b.level, d)).sum();
                    r = r.max(s);
                }
                r
            })
            .reduce(|| 0.0, f64::max);
        let mut dt = self.cfl / rate;
        if !self.sys.nc_vars().is_empty() {
            // stiffness of the non-conservative dilatation term (multiphase)
            let src = self
                .blocks
                .par_iter()
                .filter(|b| b.leaf)
                .map(|b| {
                    let mut r: f64 = 0.0;
                    let bd = [self.bd(0), self.bd(1), self.bd(2)];
                    for i in 0..bd[0] {
                        for j in 0..bd[1] {
                            for k in 0..bd[2] {
                                let c = [i as i64, j as i64, k as i64];
                                let w = self.sys.to_prim(&b.u[self.cidx([i, j, k])]);
                                let mut divu = 0.0;
                                for &d in &dims {
                                    let get = |o: i64| {
                                        let mut cc = c;
                                        cc[d] += o;
                                        if cc[d] >= 0 && (cc[d] as usize) < bd[d] {
                                            self.sys.to_prim(&b.u[self.cidx(cc.map(|v| v as usize))])
                                        } else {
                                            let g: [i64; 3] = std::array::from_fn(|e| b.idx[e] * bd[e] as i64 + cc[e]);
                                            self.sys.to_prim(&self.sample(b.level, g))
                                        }
                                    };
                                    let vi = self.sys.vel_index() + d;
                                    divu += (get(1)[vi] - get(-1)[vi]).abs() / (2.0 * self.dx(b.level, d));
                                }
                                r = r.max(self.sys.dilatation_rate(&w, divu));
                            }
                        }
                    }
                    r
                })
                .reduce(|| 0.0, f64::max);
            if src > 0.0 {
                dt = dt.min(0.5 / src);
            }
        }
        dt
    }

    /// One SSP-RK3 step of size dt (global), regridding periodically.
    pub fn step(&mut self, dt: f64) {
        self.stage_dt = dt;
        for b in self.blocks.iter_mut().filter(|b| b.leaf) {
            b.u0 = b.u.clone();
        }
        let s1 = self.compute_rhs();
        self.blocks.par_iter_mut().filter(|b| b.leaf).for_each(|b| {
            for (u, r) in b.u.iter_mut().zip(b.rhs.iter()) {
                for q in 0..N {
                    u[q] += dt * r[q];
                }
            }
        });
        let s2 = self.compute_rhs();
        self.blocks.par_iter_mut().filter(|b| b.leaf).for_each(|b| {
            for ((u, r), u0) in b.u.iter_mut().zip(b.rhs.iter()).zip(b.u0.iter()) {
                for q in 0..N {
                    u[q] = 0.75 * u0[q] + 0.25 * (u[q] + dt * r[q]);
                }
            }
        });
        let s3 = self.compute_rhs();
        self.blocks.par_iter_mut().filter(|b| b.leaf).for_each(|b| {
            for ((u, r), u0) in b.u.iter_mut().zip(b.rhs.iter()).zip(b.u0.iter()) {
                for q in 0..N {
                    u[q] = u0[q] / 3.0 + 2.0 / 3.0 * (u[q] + dt * r[q]);
                }
            }
        });
        for s in [s1, s2, s3] {
            self.stats.add(&s);
        }
        if let Some(rx) = &self.relax {
            let sys = &self.sys;
            self.blocks.par_iter_mut().filter(|b| b.leaf).for_each(|b| {
                for u in b.u.iter_mut() {
                    rx(sys, u, dt);
                }
            });
        }
        self.advance_clock(dt);
        self.regrid_if_due();
    }

    /// Advance time and the step counter (after a step's updates).
    pub fn advance_clock(&mut self, dt: f64) {
        self.t += dt;
        self.steps += 1;
    }

    /// Whether the next `regrid_if_due` will regrid.
    pub fn regrid_is_due(&self) -> bool {
        self.steps >= self.next_regrid
    }

    /// Regrid if the (adaptive) interval has elapsed; Some(changed) if it ran.
    pub fn regrid_if_due(&mut self) -> Option<bool> {
        if !self.regrid_is_due() {
            return None;
        }
        let t0 = std::time::Instant::now();
        let before = self.structure();
        self.regrid();
        let changed = self.structure() != before;
        self.schedule_regrid(changed);
        self.prof[4] += t0.elapsed().as_secs_f64();
        Some(changed)
    }

    /// (level, coordinates, leaf) of every block, in id order.
    pub fn structure(&self) -> Vec<(usize, [i64; 3], bool)> {
        self.blocks.iter().map(|b| (b.level, b.idx, b.leaf)).collect()
    }

    /// Set the next regrid step after a regrid that `changed` the grid or not.
    pub fn schedule_regrid(&mut self, changed: bool) {
        self.regrid_interval = if changed { self.cfg.regrid_every } else { (2 * self.regrid_interval).min(2 * self.cfg.regrid_every) };
        self.next_regrid = self.steps + self.regrid_interval;
    }

    /// Hierarchy description for an accelerator backend (block list, dense
    /// per-level index, ancestors to restrict, refluxing faces, operators).
    pub fn layout(&self) -> AmrLayout {
        let dims = self.dims();
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let levels = self.cfg.levels;
        let edge: Vec<Vec<usize>> = (0..levels)
            .map(|l| {
                (0..self.blocks.len())
                    .filter(|&i| {
                        let b = &self.blocks[i];
                        !b.leaf && b.level == l && self.neighbours(l, b.idx).iter().any(|nk| self.map.get(&(l as u8, *nk)).map_or(true, |&j| self.blocks[j].leaf))
                    })
                    .collect()
            })
            .collect();
        let mut reflux = Vec::new();
        for (id, b) in self.blocks.iter().enumerate() {
            let l = b.level;
            if !b.leaf || l + 1 >= levels {
                continue;
            }
            for &d in &dims {
                let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
                let (e0, e1) = (others[0], others[1]);
                let nbk = self.cfg.base_blocks[d] as i64 * (1 << l);
                for side in 0..2 {
                    let mut nk = b.idx;
                    nk[d] += if side == 1 { 1 } else { -1 };
                    if nk[d] < 0 || nk[d] >= nbk {
                        if self.bc[d][side] != Bc::Periodic {
                            continue;
                        }
                        nk[d] = nk[d].rem_euclid(nbk);
                    }
                    if !self.map.get(&(l as u8, nk)).is_some_and(|&j| !self.blocks[j].leaf) {
                        continue;
                    }
                    for a in 0..bd[e0] {
                        for bb in 0..bd[e1] {
                            let mut fine = Vec::new();
                            for fa in if dims.contains(&e0) { 0..2 } else { 0..1 } {
                                for fb in if dims.contains(&e1) { 0..2 } else { 0..1 } {
                                    let ga = 2 * (b.idx[e0] * bd[e0] as i64 + a as i64) + fa;
                                    let gb = 2 * (b.idx[e1] * bd[e1] as i64 + bb as i64) + fb;
                                    let mut ck = nk;
                                    ck[d] = 2 * nk[d] + if side == 1 { 0 } else { 1 };
                                    ck[e0] = if dims.contains(&e0) { ga.div_euclid(bd[e0] as i64) } else { 0 };
                                    ck[e1] = if dims.contains(&e1) { gb.div_euclid(bd[e1] as i64) } else { 0 };
                                    let Some(&cid) = self.map.get(&((l + 1) as u8, ck)) else { continue };
                                    if !self.blocks[cid].leaf {
                                        continue;
                                    }
                                    let la = if dims.contains(&e0) { ga.rem_euclid(bd[e0] as i64) as usize } else { 0 };
                                    let lb = if dims.contains(&e1) { gb.rem_euclid(bd[e1] as i64) as usize } else { 0 };
                                    fine.push((cid, la * bd[e1] + lb));
                                }
                            }
                            if fine.is_empty() {
                                continue;
                            }
                            let mut c = [0usize; 3];
                            c[d] = if side == 1 { bd[d] - 1 } else { 0 };
                            c[e0] = a;
                            c[e1] = bb;
                            reflux.push(RefluxFace { block: id, d, side, face: a * bd[e1] + bb, cell: self.cidx(c), fine });
                        }
                    }
                }
            }
        }
        AmrLayout {
            ng: self.ng,
            bd,
            active: std::array::from_fn(|d| self.active(d)),
            levels,
            nb0: self.map.nb0,
            level: self.blocks.iter().map(|b| b.level).collect(),
            idx: self.blocks.iter().map(|b| b.idx).collect(),
            leaf: self.blocks.iter().map(|b| b.leaf).collect(),
            map: self.map.levels.clone(),
            edge,
            reflux,
            wp: self.wp,
            wr: self.wr,
            dx: (0..levels).map(|l| std::array::from_fn(|d| self.dx(l, d))).collect(),
        }
    }

    pub fn run_to(&mut self, t_end: f64) -> usize {
        let mut n = 0;
        while self.t < t_end * (1.0 - 1e-14) {
            let dt = if self.subcycle { self.dt_cfl_subcycled() } else { self.dt_cfl() }.min(t_end - self.t);
            assert!(dt.is_finite() && dt > 0.0, "time step collapsed");
            if self.subcycle {
                self.step_subcycled(dt);
            } else {
                self.step(dt);
            }
            n += 1;
        }
        n
    }

    // ---------------------------------------------------------------- subcycling

    /// CFL time step of the blocks at level l (leaves and ancestors).
    fn dt_cfl_level(&self, l: usize) -> f64 {
        let dims = self.dims();
        let rate = self
            .blocks
            .par_iter()
            .filter(|b| b.level == l)
            .map(|b| {
                let mut r: f64 = 0.0;
                for u in &b.u {
                    let w = self.sys.to_prim(u);
                    let sp: f64 = dims.iter().map(|&d| self.sys.max_speed(&w, d) / self.dx(l, d)).sum();
                    r = r.max(sp);
                }
                r
            })
            .reduce(|| 0.0, f64::max);
        if rate > 0.0 {
            self.cfl / rate
        } else {
            f64::INFINITY
        }
    }

    /// Level-0 time step for subcycling: every level must satisfy its own
    /// CFL condition with Δt_l = Δt₀/2^l.
    pub fn dt_cfl_subcycled(&self) -> f64 {
        let lmin = self.coarsest_leaf_level();
        let mut dt = f64::INFINITY;
        for l in lmin..self.cfg.levels {
            dt = dt.min(self.dt_cfl_level(l) * (1u64 << (l - lmin)) as f64);
        }
        // the multiphase dilatation limit is enforced on the finest scale
        let g = self.dt_cfl();
        let finest = (0..self.cfg.levels).rev().find(|&l| self.blocks.iter().any(|b| b.level == l && b.leaf)).unwrap_or(0);
        dt.min(g * (1u64 << (finest - lmin)) as f64)
    }

    /// Coarsest level holding leaves: subcycling starts there (levels below
    /// are all covered and are updated by restriction only; advancing them
    /// would cost work and make the coarse step — and the regrid interval —
    /// 2^(levels) finest steps long).
    pub fn coarsest_leaf_level(&self) -> usize {
        self.blocks.iter().filter(|b| b.leaf).map(|b| b.level).min().unwrap_or(0)
    }

    /// Replace the data of all level-l blocks by their linear time
    /// interpolation at time t (saving the new state), for ghost cells of
    /// level l + 1.
    fn set_time_view(&mut self, l: usize, t: f64) {
        let (t0, dt) = (self.level_t0[l], self.level_dt[l]);
        let th = if dt > 0.0 { ((t - t0) / dt).clamp(0.0, 1.0) } else { 1.0 };
        self.blocks.par_iter_mut().filter(|b| b.level == l && !b.u0.is_empty()).for_each(|b| {
            b.u_save.clone_from(&b.u);
            for (u, u0) in b.u.iter_mut().zip(b.u0.iter()) {
                for q in 0..N {
                    u[q] = (1.0 - th) * u0[q] + th * u[q];
                }
            }
        });
    }

    fn restore_time_view(&mut self, l: usize) {
        self.blocks.par_iter_mut().filter(|b| b.level == l && !b.u_save.is_empty()).for_each(|b| {
            std::mem::swap(&mut b.u, &mut b.u_save);
        });
    }

    /// Faces of block `id` on a coarse–fine boundary: (d, side, kind) with
    /// kind 0 = neighbour refined (this is the coarse side), 1 = neighbour
    /// coarser (this is the fine side; returns the coarse neighbour block).
    fn cf_faces(&self, id: usize) -> Vec<(usize, usize, u8, Option<usize>)> {
        let b = &self.blocks[id];
        let l = b.level;
        let dims = self.dims();
        let mut out = Vec::new();
        for &d in &dims {
            let nbk = self.cfg.base_blocks[d] as i64 * (1 << l);
            for side in 0..2 {
                let mut nk = b.idx;
                nk[d] += if side == 1 { 1 } else { -1 };
                if nk[d] < 0 || nk[d] >= nbk {
                    if self.bc[d][side] != Bc::Periodic {
                        continue;
                    }
                    nk[d] = nk[d].rem_euclid(nbk);
                }
                match self.map.get(&(l as u8, nk)) {
                    Some(&j) if !self.blocks[j].leaf => out.push((d, side, 0, None)),
                    Some(_) => {}
                    None if l > 0 => {
                        let pk: [i64; 3] = std::array::from_fn(|e| if dims.contains(&e) { nk[e] >> 1 } else { 0 });
                        if let Some(&cj) = self.map.get(&((l - 1) as u8, pk)) {
                            out.push((d, side, 1, Some(cj)));
                        }
                    }
                    None => {}
                }
            }
        }
        out
    }

    /// Accumulate the stage's boundary fluxes into the registers.
    fn accumulate_registers(&mut self, ids: &[usize], wdt: f64) {
        let dims = self.dims();
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let nsub = (1usize << (dims.len() - 1)) as f64;
        for &id in ids {
            if !self.blocks[id].leaf {
                continue;
            }
            for (d, side, kind, coarse) in self.cf_faces(id) {
                let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
                let (e0, e1) = (others[0], others[1]);
                let face = self.blocks[id].face[d][side].clone();
                if kind == 0 {
                    for (fi, f) in face.iter().enumerate() {
                        let e = self.creg.entry((id, d as u8, side as u8, fi)).or_insert([0.0; N]);
                        for q in 0..N {
                            e[q] += wdt * f[q];
                        }
                    }
                } else {
                    // fine face (a, b) → coarse face index on the coarse block's opposite side
                    let cj = coarse.unwrap();
                    let b = &self.blocks[id];
                    let cb = &self.blocks[cj];
                    for a in 0..bd[e0] {
                        for bb in 0..bd[e1] {
                            let ga = b.idx[e0] * bd[e0] as i64 + a as i64;
                            let gb = b.idx[e1] * bd[e1] as i64 + bb as i64;
                            let ca = if dims.contains(&e0) { (ga.div_euclid(2) - cb.idx[e0] * bd[e0] as i64) as usize } else { 0 };
                            let cbb = if dims.contains(&e1) { (gb.div_euclid(2) - cb.idx[e1] * bd[e1] as i64) as usize } else { 0 };
                            if ca >= bd[e0] || cbb >= bd[e1] {
                                continue;
                            }
                            let f = face[a * bd[e1] + bb];
                            let key = (cj, d as u8, (1 - side) as u8, ca * bd[e1] + cbb);
                            let e = self.freg.entry(key).or_insert([0.0; N]);
                            for q in 0..N {
                                e[q] += wdt * f[q] / nsub;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Correct coarse leaves of level l with the flux registers.
    fn reflux_level(&mut self, l: usize) {
        let dims = self.dims();
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let keys: Vec<(usize, u8, u8, usize)> = self.creg.keys().filter(|k| self.blocks[k.0].level == l).cloned().collect();
        for key in keys {
            let (id, d, side, fi) = key;
            let Some(fr) = self.freg.get(&key).copied() else { continue };
            let cr = self.creg[&key];
            let d = d as usize;
            let others: Vec<usize> = (0..3).filter(|&e| e != d).collect();
            let (_e0, e1) = (others[0], others[1]);
            let (a, bb) = (fi / bd[e1], fi % bd[e1]);
            let mut c = [0usize; 3];
            c[d] = if side == 1 { bd[d] - 1 } else { 0 };
            c[others[0]] = a;
            c[e1] = bb;
            let dx = self.dx(l, d);
            let cell = self.cidx(c);
            for q in 0..N {
                let df = fr[q] - cr[q];
                self.blocks[id].u[cell][q] += if side == 1 { -df / dx } else { df / dx };
            }
        }
        let _ = dims;
    }

    /// Advance all blocks of level l (leaves and covered ancestors) by dt,
    /// then the finer level twice by dt/2, then reflux and restrict.
    fn advance_level(&mut self, l: usize, dt: f64) {
        let ids: Vec<usize> = (0..self.blocks.len()).filter(|&i| self.blocks[i].level == l).collect();
        if ids.is_empty() {
            return;
        }
        let t0 = self.level_t0[l];
        self.level_dt[l] = dt;
        for &id in &ids {
            let b = &mut self.blocks[id];
            b.u0 = b.u.clone();
        }
        // registers of this level's coarse faces / fine contributions to them
        let lvl_ids: std::collections::HashSet<usize> = ids.iter().cloned().collect();
        self.creg.retain(|k, _| !lvl_ids.contains(&k.0));
        self.freg.retain(|k, _| !lvl_ids.contains(&k.0));
        self.stage_dt = dt;
        let stage_t = [t0, t0 + dt, t0 + 0.5 * dt];
        let weights = [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0];
        let ghost_ns = std::sync::atomic::AtomicU64::new(0);
        let flux_ns = std::sync::atomic::AtomicU64::new(0);
        for s in 0..3 {
            if l > 0 {
                self.set_time_view(l - 1, stage_t[s]);
            }
            let res: Vec<(Vec<V<N>>, [[Vec<V<N>>; 2]; 3], Stats)> = ids.par_iter().map(|&id| self.block_rhs(id, &ghost_ns, &flux_ns)).collect();
            if l > 0 {
                self.restore_time_view(l - 1);
            }
            for (&id, (rhs, face, st)) in ids.iter().zip(res) {
                self.blocks[id].rhs = rhs;
                self.blocks[id].face = face;
                self.stats.add(&st);
            }
            self.accumulate_registers(&ids, weights[s] * dt);
            let (a0, a1) = [(0.0, 1.0), (0.75, 0.25), (1.0 / 3.0, 2.0 / 3.0)][s];
            for &id in &ids {
                let b = &mut self.blocks[id];
                let (u, r, u0) = (&mut b.u, &b.rhs, &b.u0);
                for ((uc, rc), u0c) in u.iter_mut().zip(r.iter()).zip(u0.iter()) {
                    for q in 0..N {
                        uc[q] = a0 * u0c[q] + a1 * (uc[q] + dt * rc[q]);
                    }
                }
            }
        }
        let nthreads = rayon::current_num_threads() as f64;
        self.prof[1] += ghost_ns.into_inner() as f64 * 1e-9 / nthreads;
        self.prof[2] += flux_ns.into_inner() as f64 * 1e-9 / nthreads;
        if l + 1 < self.cfg.levels && self.blocks.iter().any(|b| b.level == l + 1) {
            self.level_t0[l + 1] = t0;
            self.advance_level(l + 1, 0.5 * dt);
            self.level_t0[l + 1] = t0 + 0.5 * dt;
            self.advance_level(l + 1, 0.5 * dt);
            let tr = std::time::Instant::now();
            self.reflux_level(l);
            self.prof[3] += tr.elapsed().as_secs_f64();
            // covered cells take the restriction of the fine solution
            let tq = std::time::Instant::now();
            self.restrict_level_full(l);
            self.prof[0] += tq.elapsed().as_secs_f64();
        }
        if let Some(rx) = &self.relax {
            let sys = &self.sys;
            self.blocks.par_iter_mut().filter(|b| b.leaf && b.level == l).for_each(|b| {
                for u in b.u.iter_mut() {
                    rx(sys, u, dt);
                }
            });
        }
        self.level_t0[l] = t0 + dt;
    }

    /// Restrict every ancestor block of level l from level l + 1.
    fn restrict_level_full(&mut self, l: usize) {
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        let ids: Vec<usize> = (0..self.blocks.len()).filter(|&i| !self.blocks[i].leaf && self.blocks[i].level == l).collect();
        let new: Vec<Vec<V<N>>> = ids.par_iter().map(|&id| self.restrict_box(l, self.blocks[id].idx, [0; 3], bd)).collect();
        for (id, u) in ids.into_iter().zip(new) {
            self.blocks[id].u = u;
        }
    }

    /// One subcycled step: level 0 by dt (finer levels by dt/2^l).
    pub fn step_subcycled(&mut self, dt: f64) {
        // registers are keyed by block id: regridding reorders blocks
        self.creg.clear();
        self.freg.clear();
        for l in 0..self.cfg.levels {
            self.level_t0[l] = self.t;
        }
        // ancestors must hold consistent data before they are advanced
        for l in (0..self.cfg.levels.saturating_sub(1)).rev() {
            self.restrict_level_full(l);
        }
        let lmin = self.coarsest_leaf_level();
        self.advance_level(lmin, dt);
        for l in (0..lmin).rev() {
            self.restrict_level_full(l);
        }
        self.t += dt;
        self.steps += 1;
        // one coarse step is 2^(L−1) finest steps: regrid every coarse step
        // (twice that at most when nothing changes) so features cannot
        // outrun the refinement buffer
        if self.next_regrid > self.steps + 2 {
            self.next_regrid = self.steps + 1;
        }
        if self.steps >= self.next_regrid {
            let t0 = std::time::Instant::now();
            let before: Vec<(usize, [i64; 3], bool)> = self.blocks.iter().map(|b| (b.level, b.idx, b.leaf)).collect();
            self.regrid();
            let changed = before.len() != self.blocks.len() || self.blocks.iter().zip(before.iter()).any(|(b, o)| (b.level, b.idx, b.leaf) != *o);
            self.regrid_interval = if changed { self.cfg.regrid_every } else { (2 * self.regrid_interval).min(2 * self.cfg.regrid_every) };
            self.next_regrid = self.steps + if changed { 1 } else { 2 };
            self.prof[4] += t0.elapsed().as_secs_f64();
        }
    }


    /// Number of leaf cells.
    pub fn cells(&self) -> usize {
        self.blocks.iter().filter(|b| b.leaf).count() * self.bcells()
    }

    /// Cells of a uniform grid at the finest level.
    pub fn uniform_cells(&self) -> usize {
        let l = self.cfg.levels - 1;
        (0..3).map(|d| self.n_level(l, d) as usize).product()
    }

    /// Leaf blocks per level.
    pub fn level_counts(&self) -> Vec<usize> {
        let mut c = vec![0; self.cfg.levels];
        for b in self.blocks.iter().filter(|b| b.leaf) {
            c[b.level] += 1;
        }
        c
    }

    /// Visit every leaf cell: (level, centre, conserved state, cell volume).
    pub fn for_each_leaf_cell(&self, mut f: impl FnMut(usize, [f64; 3], &V<N>, f64)) {
        let dims = self.dims();
        let bd = [self.bd(0), self.bd(1), self.bd(2)];
        for b in self.blocks.iter().filter(|b| b.leaf) {
            let vol: f64 = dims.iter().map(|&d| self.dx(b.level, d)).product();
            for i in 0..bd[0] {
                for j in 0..bd[1] {
                    for k in 0..bd[2] {
                        let g = [b.idx[0] * bd[0] as i64 + i as i64, b.idx[1] * bd[1] as i64 + j as i64, b.idx[2] * bd[2] as i64 + k as i64];
                        f(b.level, self.x(b.level, g), &b.u[self.cidx([i, j, k])], vol);
                    }
                }
            }
        }
    }

    /// Totals of conserved variables (Σ u ΔV over leaves).
    pub fn totals(&self) -> V<N> {
        let dims = self.dims();
        let mut s = [0.0; N];
        for b in self.blocks.iter().filter(|b| b.leaf) {
            let vol: f64 = dims.iter().map(|&d| self.dx(b.level, d)).product();
            for u in &b.u {
                for q in 0..N {
                    s[q] += u[q] * vol;
                }
            }
        }
        s
    }
}
