//! Banded LU with partial pivoting, for the Newton systems of the device
//! solver (unknowns ordered so that couplings stay within a band).
//!
//! Row r stores columns [r − kl − NB, r + kl + ku]: the extra kl columns
//! on the right hold the fill of partial pivoting (U has upper bandwidth
//! kl + ku), the NB on the left the multipliers of rows swapped within a
//! panel. Rows are scaled to unit maximum before factoring (the drift–
//! diffusion rows differ by many orders of magnitude).
//!
//! Blocked right-looking factorization (as LAPACK's dgbtrf): a panel of NB
//! columns is factored with updates restricted to the panel, then its
//! effect on the trailing rows is applied at once as a rank-NB update —
//! long contiguous vector loops, parallel over rows. Pivoting is identical
//! to the unblocked algorithm.

const NB: usize = 32;

/// dst −= l·src, compiled for AVX2 when the CPU has it (same rounding as
/// the scalar loop: no fused multiply–add).
#[inline]
fn axpy(dst: &mut [f64], l: f64, src: &[f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU supports AVX2
            unsafe { axpy_avx2(dst, l, src) };
            return;
        }
    }
    for (x, u) in dst.iter_mut().zip(src) {
        *x -= l * u;
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn axpy_avx2(dst: &mut [f64], l: f64, src: &[f64]) {
    for (x, u) in dst.iter_mut().zip(src) {
        *x -= l * u;
    }
}

use rayon::prelude::*;

pub struct Band {
    pub n: usize,
    pub kl: usize,
    pub ku: usize,
    w: usize,
    a: Vec<f64>,
}

impl Band {
    pub fn new(n: usize, kl: usize, ku: usize) -> Self {
        let w = 2 * kl + ku + NB + 1;
        Band { n, kl, ku, w, a: vec![0.0; n * w] }
    }
    pub fn clear(&mut self) {
        self.a.iter_mut().for_each(|x| *x = 0.0);
    }
    #[inline]
    fn at(&self, r: usize, c: usize) -> usize {
        debug_assert!(c + self.kl + NB >= r && c <= r + self.kl + self.ku);
        r * self.w + (c + self.kl + NB - r)
    }
    /// A[r][c] += v (|c − r| within the original band).
    #[inline]
    pub fn add(&mut self, r: usize, c: usize, v: f64) {
        assert!(c + self.kl >= r && c <= r + self.ku, "entry ({r},{c}) outside the band");
        let i = self.at(r, c);
        self.a[i] += v;
    }
    pub fn get(&self, r: usize, c: usize) -> f64 {
        if c + self.kl + NB < r || c > r + self.kl + self.ku {
            0.0
        } else {
            self.a[self.at(r, c)]
        }
    }

    /// y = A x (before factoring).
    pub fn mul(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n)
            .map(|r| {
                let c0 = r.saturating_sub(self.kl);
                let c1 = (r + self.ku).min(self.n - 1);
                (c0..=c1).map(|c| self.get(r, c) * x[c]).sum()
            })
            .collect()
    }

    /// Solve A x = b in place (A is destroyed). Returns None if singular.
    pub fn solve(&mut self, b: &mut [f64]) -> Option<()> {
        let (n, kl, ku, w) = (self.n, self.kl, self.ku, self.w);
        // row equilibration
        for r in 0..n {
            let row = &mut self.a[r * w..(r + 1) * w];
            let m = row.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            if m == 0.0 || !m.is_finite() {
                return None;
            }
            row.iter_mut().for_each(|v| *v /= m);
            b[r] /= m;
        }
        let prof = std::env::var("BAND_PROF").is_ok();
        let t0 = std::time::Instant::now();
        let (mut tp, mut tu, mut tt) = (0.0, 0.0, 0.0);
        let off = kl + NB; // column c of row r at r·w + c + off − r
        let mut c0 = 0;
        while c0 < n {
            let c1 = (c0 + NB).min(n);
            let ta = std::time::Instant::now();
            // panel: pivot, swap, eliminate within the panel columns
            for c in c0..c1 {
                let rmax = (c + kl).min(n - 1);
                let mut p = c;
                let mut best = self.a[c * w + off].abs();
                for r in c + 1..=rmax {
                    let v = self.a[r * w + c + off - r].abs();
                    if v > best {
                        best = v;
                        p = r;
                    }
                }
                if best == 0.0 || !best.is_finite() {
                    return None;
                }
                let cend = (c + kl + ku).min(n - 1);
                if p != c {
                    let (lo, hi) = self.a.split_at_mut(p * w);
                    lo[c * w + c0 + off - c..=c * w + cend + off - c].swap_with_slice(&mut hi[c0 + off - p..=cend + off - p]);
                    b.swap(c, p);
                }
                let d = self.a[c * w + off];
                // row c has no entries beyond c + kl + ku
                let jend = c1.min(cend + 1);
                let (lo, hi) = self.a.split_at_mut((c + 1) * w);
                let prow = &lo[c * w + c + 1 + off - c..c * w + jend + off - c];
                for r in c + 1..=rmax {
                    let row = &mut hi[(r - c - 1) * w..(r - c) * w];
                    let l = row[c + off - r] / d;
                    if l == 0.0 {
                        continue;
                    }
                    row[c + off - r] = l;
                    axpy(&mut row[c + 1 + off - r..jend + off - r], l, prow);
                    b[r] -= l * b[c];
                }
            }
            tp += ta.elapsed().as_secs_f64();
            if c1 == n {
                break;
            }
            let ta = std::time::Instant::now();
            // U12: apply the panel's eliminations to the trailing columns of
            // the panel rows — a triangular solve, independent per column:
            // parallel over column chunks
            let jmax = (c1 - 1 + kl + ku).min(n - 1);
            if jmax >= c1 {
                struct P(*mut f64);
                unsafe impl Send for P {}
                unsafe impl Sync for P {}
                let base = P(self.a.as_mut_ptr());
                let chunk = 32;
                let nchunks = (jmax + 1 - c1).div_ceil(chunk);
                let task = |q: usize| {
                    let base = &base;
                    let j0 = c1 + q * chunk;
                    let j1 = (j0 + chunk).min(jmax + 1);
                    for k in c0..c1 {
                        for k2 in c0..k {
                            // SAFETY: tasks write disjoint column ranges of
                            // rows k and only read rows k2 < k in the same
                            // columns, which no other task writes
                            unsafe {
                                let l = *base.0.add(k * w + k2 + off - k);
                                let hi = j1.min(k2 + kl + ku + 1);
                                if l == 0.0 || hi <= j0 {
                                    continue;
                                }
                                let dst = std::slice::from_raw_parts_mut(base.0.add(k * w + j0 + off - k), hi - j0);
                                let src = std::slice::from_raw_parts(base.0.add(k2 * w + j0 + off - k2), hi - j0);
                                axpy(dst, l, src);
                            }
                        }
                    }
                };
                if nchunks > 1 && NB * NB * (jmax + 1 - c1) > 100_000 {
                    (0..nchunks).into_par_iter().for_each(task);
                } else {
                    (0..nchunks).for_each(task);
                }
            }
            tu += ta.elapsed().as_secs_f64();
            let ta = std::time::Instant::now();
            // trailing rows: A22 −= L21 U12 (parallel over rows)
            let rend = (c1 - 1 + kl).min(n - 1);
            if rend >= c1 {
                let (top, bottom) = self.a.split_at_mut(c1 * w);
                let top = &*top;
                let rows = &mut bottom[..(rend + 1 - c1) * w];
                let update = |(ri, row): (usize, &mut [f64])| {
                    let r = c1 + ri;
                    for k in c0..c1 {
                        if k + kl + NB < r {
                            continue;
                        }
                        let l = row[k + off - r];
                        if l == 0.0 {
                            continue;
                        }
                        let jend = (k + kl + ku).min(n - 1);
                        if jend < c1 {
                            continue;
                        }
                        let urow = &top[k * w + c1 + off - k..=k * w + jend + off - k];
                        axpy(&mut row[c1 + off - r..=jend + off - r], l, urow);
                    }
                };
                if (rend + 1 - c1) * NB * (kl + ku) > 50_000 {
                    rows.par_chunks_mut(w).enumerate().for_each(update);
                } else {
                    rows.chunks_mut(w).enumerate().for_each(update);
                }
            }
            tt += ta.elapsed().as_secs_f64();
            c0 = c1;
        }
        if prof {
            eprintln!("band n {n} kl {kl}: total {:.3} s (panel {tp:.3}, U12 {tu:.3}, trailing {tt:.3}, equilibration+rest {:.3})", t0.elapsed().as_secs_f64(), t0.elapsed().as_secs_f64() - tp - tu - tt);
        }
        // back substitution with U
        for r in (0..n).rev() {
            let cend = (r + kl + ku).min(n - 1);
            let mut s = b[r];
            for c in r + 1..=cend {
                s -= self.a[r * w + c + off - r] * b[c];
            }
            b[r] = s / self.a[r * w + off];
        }
        Some(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn random_band() {
        for (n, kl, ku) in [(300, 7, 5), (1000, 40, 40), (517, 100, 3), (64, 63, 63), (20, 2, 9)] {
            check(n, kl, ku);
        }
    }

    fn check(n: usize, kl: usize, ku: usize) {
        let mut s = 12345u64;
        let mut rnd = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s % 10000) as f64 / 5000.0 - 1.0
        };
        let mut a = Band::new(n, kl, ku);
        for r in 0..n {
            for c in r.saturating_sub(kl)..=(r + ku).min(n - 1) {
                // weak diagonal: pivoting needed
                a.add(r, c, rnd() * if r == c { 0.01 } else { 1.0 });
            }
        }
        let x: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
        let b0 = a.mul(&x);
        let mut b = b0.clone();
        let keep: Vec<(usize, usize, f64)> = (0..n).flat_map(|r| (r.saturating_sub(kl)..=(r + ku).min(n - 1)).map(move |c| (r, c))).map(|(r, c)| (r, c, a.get(r, c))).collect();
        a.solve(&mut b).unwrap();
        // backward error (the random matrices can be ill-conditioned)
        let mut res = b0.clone();
        for &(r, c, v) in &keep {
            res[r] -= v * b[c];
        }
        let scale = keep.iter().map(|&(_, c, v)| (v * b[c]).abs()).fold(0.0, f64::max);
        let err = res.iter().fold(0.0f64, |m, v| m.max(v.abs())) / scale;
        assert!(err < 1e-12, "n {n} kl {kl} ku {ku}: backward error {err}");
    }
}
