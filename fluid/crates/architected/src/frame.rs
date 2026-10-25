//! A small direct 2D frame solver for benchmark beams.
//!
//! Used for the exact cantilever and fixed-fixed checks that calibrate the
//! element implementation before it is trusted inside a unit cell.

use crate::beam::{frame2d_global, frame2d_local, timoshenko_phi};
use crate::linalg::{solve, Dense};

/// 2D frame member.
#[derive(Clone, Copy, Debug)]
pub struct Member {
    pub a: usize,
    pub b: usize,
    pub e: f64,
    pub area: f64,
    pub inertia: f64,
}

/// 2D frame with supports and a linear-elastic solve.
#[derive(Clone, Debug)]
pub struct Frame {
    pub nodes: Vec<[f64; 2]>,
    pub members: Vec<Member>,
    /// `[u_x, u_y, theta_z]` fixed flags per node.
    pub supports: Vec<[bool; 3]>,
    pub nu: Option<f64>,
}

impl Frame {
    pub fn new(nodes: Vec<[f64; 2]>, members: Vec<Member>) -> Self {
        let n = nodes.len();
        Self {
            nodes,
            members,
            supports: vec![[false; 3]; n],
            nu: None,
        }
    }

    /// A straight chain of `n` members along `x` from `0` to `length`.
    pub fn straight_beam(length: f64, n: usize, e: f64, area: f64, inertia: f64) -> Self {
        let nodes: Vec<[f64; 2]> = (0..=n)
            .map(|i| [length * i as f64 / n as f64, 0.0])
            .collect();
        let members: Vec<Member> = (0..n)
            .map(|i| Member {
                a: i,
                b: i + 1,
                e,
                area,
                inertia,
            })
            .collect();
        Self::new(nodes, members)
    }

    pub fn assemble(&self) -> Dense {
        let n = 3 * self.nodes.len();
        let mut k = Dense::zeros(n);
        for m in &self.members {
            let pa = self.nodes[m.a];
            let pb = self.nodes[m.b];
            let dx = pb[0] - pa[0];
            let dy = pb[1] - pa[1];
            let l = (dx * dx + dy * dy).sqrt();
            let (c, s) = (dx / l, dy / l);
            let phi = match self.nu {
                Some(nu) => timoshenko_phi(l, m.e, m.area, m.inertia, nu),
                None => 0.0,
            };
            let kl = frame2d_local(l, m.e, m.area, m.inertia, phi);
            let kg = frame2d_global(&kl, c, s);
            let idx = [
                3 * m.a,
                3 * m.a + 1,
                3 * m.a + 2,
                3 * m.b,
                3 * m.b + 1,
                3 * m.b + 2,
            ];
            for (ii, &ri) in idx.iter().enumerate() {
                for (jj, &rj) in idx.iter().enumerate() {
                    k.add(ri, rj, kg[ii][jj]);
                }
            }
        }
        k
    }

    /// Solve `K u = f`. Returns the nodal displacements and the reaction
    /// forces at the supports.
    pub fn solve(&self, loads: &[[f64; 3]]) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let n = 3 * self.nodes.len();
        assert_eq!(loads.len(), self.nodes.len());
        let mut full_f = vec![0.0; n];
        for (i, l) in loads.iter().enumerate() {
            full_f[3 * i] = l[0];
            full_f[3 * i + 1] = l[1];
            full_f[3 * i + 2] = l[2];
        }
        let fixed: Vec<bool> = self
            .supports
            .iter()
            .flat_map(|s| s.iter().copied())
            .collect();
        let free_idx: Vec<usize> = (0..n).filter(|&i| !fixed[i]).collect();
        let mut map = vec![usize::MAX; n];
        for (new, &old) in free_idx.iter().enumerate() {
            map[old] = new;
        }
        let kfull = self.assemble();
        let nr = free_idx.len();
        let mut kr = Dense::zeros(nr);
        for (i, &oi) in free_idx.iter().enumerate() {
            for (j, &oj) in free_idx.iter().enumerate() {
                kr.set(i, j, kfull.get(oi, oj));
            }
        }
        let mut fr: Vec<f64> = free_idx.iter().map(|&i| full_f[i]).collect();
        solve(&mut kr, &mut fr).expect("singular frame");
        let mut u = vec![0.0; n];
        for (i, &oi) in free_idx.iter().enumerate() {
            u[oi] = fr[i];
        }
        // Reactions: r = K u - f at the fixed DOFs.
        let ku = kfull.mul_vec(&u);
        let mut reactions = vec![[0.0; 3]; self.nodes.len()];
        for i in 0..n {
            if fixed[i] {
                let r = ku[i] - full_f[i];
                reactions[i / 3][i % 3] = r;
            }
        }
        let mut disp = vec![[0.0; 3]; self.nodes.len()];
        for i in 0..self.nodes.len() {
            disp[i] = [u[3 * i], u[3 * i + 1], u[3 * i + 2]];
        }
        (disp, reactions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cantilever_tip_deflection_is_exact() {
        // Exact Euler-Bernoulli tip load: delta = P L^3 / (3 E I).
        let (e, b, h, l, p) = (200e9, 0.02, 0.01, 1.0, 100.0);
        let area = b * h;
        let inertia = b * h * h * h / 12.0;
        // One element is enough for exact nodal displacement under end load.
        let mut f = Frame::straight_beam(l, 1, e, area, inertia);
        f.supports[0] = [true, true, true];
        let mut loads = vec![[0.0; 3]; f.nodes.len()];
        let tip = f.nodes.len() - 1;
        loads[tip][1] = -p;
        let (u, _) = f.solve(&loads);
        let exact = p * l * l * l / (3.0 * e * inertia);
        assert!((u[tip][1].abs() - exact).abs() / exact < 1e-12);
    }
    #[test]
    fn timoshenko_adds_the_exact_shear_deflection() {
        // A stubby Timoshenko cantilever with an end load deflects by the
        // Euler-Bernoulli value plus P L/(G A_s). One element is exact for
        // this end-load problem.
        let (e, nu, l, p) = (70e9, 0.3, 1.0, 100.0);
        let (b, h) = (0.3, 0.3);
        let area = b * h;
        let inertia = b * h * h * h / 12.0;
        let mut f = Frame::straight_beam(l, 1, e, area, inertia);
        f.nu = Some(nu);
        f.supports[0] = [true, true, true];
        let tip = f.nodes.len() - 1;
        let mut loads = vec![[0.0; 3]; f.nodes.len()];
        loads[tip][1] = -p;
        let (u, _) = f.solve(&loads);
        let g = e / (2.0 * (1.0 + nu));
        let a_s = (5.0 / 6.0) * area;
        let exact = p * l * l * l / (3.0 * e * inertia) + p * l / (g * a_s);
        assert!((u[tip][1].abs() - exact).abs() / exact < 1e-12);
    }

    #[test]
    fn fixed_fixed_central_deflection_is_exact() {
        // Exact fixed-fixed central point load: delta = P L^3 / (192 E I).
        let (e, b, h, l, p) = (200e9, 0.02, 0.01, 1.0, 100.0);
        let area = b * h;
        let inertia = b * h * h * h / 12.0;
        let mut f = Frame::straight_beam(l, 2, e, area, inertia);
        f.supports[0] = [true, true, true];
        f.supports[2] = [true, true, true];
        let mut loads = vec![[0.0; 3]; f.nodes.len()];
        loads[1][1] = -p;
        let (u, _) = f.solve(&loads);
        let exact = p * l * l * l / (192.0 * e * inertia);
        assert!((u[1][1].abs() - exact).abs() / exact < 1e-12);
    }
}
