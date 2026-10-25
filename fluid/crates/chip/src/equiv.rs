//! Is the netlist the design? An RTL interpreter and an equivalence check.
//!
//! The interpreter executes the Verilog directly on machine words (it takes
//! only the branch that is true, adds with the CPU's adder), sharing nothing
//! with synthesis but the width rules, which are the language's. The
//! netlist is evaluated gate by gate. Both are put in the same state with
//! the same inputs, and their outputs and next states compared. Registers
//! map one-to-one, so agreement on every (input, state) pair is a proof
//! that the circuits are equivalent: when the space has at most
//! 2^`EXHAUSTIVE_BITS` points it is enumerated completely; otherwise random
//! points are checked and the result says so.

use std::collections::BTreeMap;

use crate::rtl::{Dir, Expr, Module, Stmt};
use crate::synth::{width, Netlist};

pub const EXHAUSTIVE_BITS: usize = 22;

fn mask(w: usize) -> u128 {
    if w >= 128 { u128::MAX } else { (1u128 << w) - 1 }
}

/// Values of the inputs and registers (by name) for one evaluation.
pub struct Interp<'a> {
    pub m: &'a Module,
    pub vals: BTreeMap<String, u128>,
    wires: BTreeMap<String, u128>,
}

impl<'a> Interp<'a> {
    pub fn new(m: &'a Module, vals: BTreeMap<String, u128>) -> Self {
        Interp { m, vals, wires: BTreeMap::new() }
    }

    fn signal(&mut self, n: &str) -> Result<u128, String> {
        if let Some(&v) = self.vals.get(n) {
            return Ok(v);
        }
        if let Some(&v) = self.wires.get(n) {
            return Ok(v);
        }
        let w = self.m.widths[n];
        let mut v = 0u128;
        for (name, part, e) in self.m.assigns.clone() {
            if name != n {
                continue;
            }
            let (msb, lsb) = part.unwrap_or((w - 1, 0));
            let pw = msb - lsb + 1;
            let x = self.eval(&e, pw.max(width(self.m, &e)?))? & mask(pw);
            v |= x << lsb;
        }
        self.wires.insert(n.to_string(), v);
        Ok(v)
    }

    pub fn eval(&mut self, e: &Expr, w: usize) -> Result<u128, String> {
        let m = self.m;
        let r = match e {
            Expr::Const(c, _) => *c,
            Expr::Ident(n) => self.signal(n)?,
            Expr::Slice(x, msb, lsb) => (self.eval(x, width(m, x)?)? >> lsb) & mask(msb - lsb + 1),
            Expr::Concat(parts) => {
                let mut v = 0u128;
                for p in parts {
                    let pw = width(m, p)?;
                    v = (v << pw) | (self.eval(p, pw)? & mask(pw));
                }
                v
            }
            Expr::Repeat(n, x) => {
                let xw = width(m, x)?;
                let one = self.eval(x, xw)? & mask(xw);
                (0..*n).fold(0u128, |acc, _| (acc << xw) | one)
            }
            Expr::Unary(op, a) => match *op {
                "~" => !self.eval(a, w)?,
                "-" => 0u128.wrapping_sub(self.eval(a, w)?),
                "!" => (self.eval(a, width(m, a)?)? & mask(width(m, a)?) == 0) as u128,
                "|" => (self.eval(a, width(m, a)?)? & mask(width(m, a)?) != 0) as u128,
                "&" => {
                    let aw = width(m, a)?;
                    (self.eval(a, aw)? & mask(aw) == mask(aw)) as u128
                }
                "^" => {
                    let aw = width(m, a)?;
                    ((self.eval(a, aw)? & mask(aw)).count_ones() & 1) as u128
                }
                _ => unreachable!(),
            },
            Expr::Binary(op, a, b) => match *op {
                "&" => self.eval(a, w)? & self.eval(b, w)?,
                "|" => self.eval(a, w)? | self.eval(b, w)?,
                "^" => self.eval(a, w)? ^ self.eval(b, w)?,
                "+" => self.eval(a, w)?.wrapping_add(self.eval(b, w)?),
                "-" => self.eval(a, w)?.wrapping_sub(self.eval(b, w)?),
                "*" => self.eval(a, w)?.wrapping_mul(self.eval(b, w)?),
                "==" | "!=" | "<" | "<=" | ">" | ">=" => {
                    let ow = width(m, a)?.max(width(m, b)?);
                    let (x, y) = (self.eval(a, ow)? & mask(ow), self.eval(b, ow)? & mask(ow));
                    (match *op {
                        "==" => x == y,
                        "!=" => x != y,
                        "<" => x < y,
                        "<=" => x <= y,
                        ">" => x > y,
                        _ => x >= y,
                    }) as u128
                }
                "&&" | "||" => {
                    let x = self.eval(a, width(m, a)?)? & mask(width(m, a)?) != 0;
                    let y = self.eval(b, width(m, b)?)? & mask(width(m, b)?) != 0;
                    (if *op == "&&" { x && y } else { x || y }) as u128
                }
                "<<" | ">>" => {
                    let x = self.eval(a, w)? & mask(w);
                    let s = self.eval(b, width(m, b)?)? & mask(width(m, b)?);
                    if s >= w as u128 {
                        0
                    } else if *op == "<<" {
                        x << s
                    } else {
                        x >> s
                    }
                }
                other => return Err(format!("operator `{other}`")),
            },
            Expr::Cond(c, a, b) => {
                let cw = width(m, c)?;
                if self.eval(c, cw)? & mask(cw) != 0 { self.eval(a, w)? } else { self.eval(b, w)? }
            }
        };
        Ok(r & mask(w))
    }

    pub fn exec(&mut self, stmts: &[Stmt], next: &mut BTreeMap<String, u128>) -> Result<(), String> {
        for st in stmts {
            match st {
                Stmt::Assign(n, part, e) => {
                    let total = self.m.widths[n];
                    let (msb, lsb) = part.unwrap_or((total - 1, 0));
                    let w = msb - lsb + 1;
                    let v = self.eval(e, w.max(width(self.m, e)?))? & mask(w);
                    let old = next[n];
                    next.insert(n.clone(), (old & !(mask(w) << lsb)) | (v << lsb));
                }
                Stmt::If(c, a, b) => {
                    let cw = width(self.m, c)?;
                    if self.eval(c, cw)? != 0 {
                        self.exec(a, next)?;
                    } else {
                        self.exec(b, next)?;
                    }
                }
                Stmt::Case(sel, arms, default) => {
                    let sw = width(self.m, sel)?;
                    let mut hit = None;
                    'arms: for (labels, body) in arms {
                        for l in labels {
                            let ow = sw.max(width(self.m, l)?);
                            if self.eval(sel, ow)? == self.eval(l, ow)? {
                                hit = Some(body);
                                break 'arms;
                            }
                        }
                    }
                    self.exec(hit.unwrap_or(default), next)?;
                }
            }
        }
        Ok(())
    }

    /// Outputs and next register state.
    pub fn step(&mut self) -> Result<(BTreeMap<String, u128>, BTreeMap<String, u128>), String> {
        let mut outs = BTreeMap::new();
        for p in self.m.ports.iter().filter(|p| p.dir == Dir::Output) {
            outs.insert(p.name.clone(), self.signal(&p.name)? & mask(p.width));
        }
        let mut next: BTreeMap<String, u128> = self.m.regs.iter().map(|r| (r.clone(), self.vals[r])).collect();
        let clocked = self.m.clocked.clone();
        self.exec(&clocked, &mut next)?;
        Ok((outs, next))
    }
}

/// Netlist outputs and next state for the same named values.
pub fn netlist_step(n: &Netlist, m: &Module, vals: &BTreeMap<String, u128>) -> (BTreeMap<String, u128>, BTreeMap<String, u128>) {
    let mut iv = vec![false; n.comb.inputs.len()];
    let pos: BTreeMap<usize, usize> = n.comb.inputs.iter().enumerate().map(|(i, &net)| (net, i)).collect();
    for (name, bits) in &n.inputs {
        for (k, net) in bits.iter().enumerate() {
            iv[pos[net]] = (vals[name] >> k) & 1 == 1;
        }
    }
    let mut fi = 0;
    for r in &m.regs {
        for k in 0..m.widths[r] {
            iv[pos[&n.flops[fi].q]] = (vals[r] >> k) & 1 == 1;
            fi += 1;
        }
    }
    iv[pos[&n.one]] = true;
    let v = n.comb.eval(&iv);
    let word = |bits: &[usize]| bits.iter().enumerate().fold(0u128, |acc, (k, &net)| acc | ((v[net] as u128) << k));
    let outs = n.outputs.iter().map(|(name, bits)| (name.clone(), word(bits))).collect();
    let mut next = BTreeMap::new();
    let mut fi = 0;
    for r in &m.regs {
        let w = m.widths[r];
        let d: Vec<usize> = n.flops[fi..fi + w].iter().map(|f| f.d).collect();
        next.insert(r.clone(), word(&d));
        fi += w;
    }
    (outs, next)
}

#[derive(Clone, Debug)]
pub struct Equivalence {
    pub points: u64,
    /// Every (input, state) point was checked: a proof.
    pub exhaustive: bool,
    pub mismatch: Option<String>,
}

/// Compare outputs and next states over the (input, register) space.
pub fn check(m: &Module, n: &Netlist, samples: u64) -> Result<Equivalence, String> {
    let mut fields: Vec<(String, usize)> = m.ports.iter().filter(|p| p.dir == Dir::Input && Some(&p.name) != m.clock.as_ref()).map(|p| (p.name.clone(), p.width)).collect();
    fields.extend(m.regs.iter().map(|r| (r.clone(), m.widths[r])));
    let total: usize = fields.iter().map(|f| f.1).sum();
    let exhaustive = total <= EXHAUSTIVE_BITS;
    let points = if exhaustive { 1u64 << total } else { samples };
    let mut rng = 0x2545F4914F6CDD1Du64;
    let mut next_rand = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    for p in 0..points {
        let mut vals = BTreeMap::new();
        let mut off = 0;
        for (name, w) in &fields {
            let v = if exhaustive {
                ((p as u128) >> off) & mask(*w)
            } else {
                (((next_rand() as u128) << 64) | next_rand() as u128) & mask(*w)
            };
            off += w;
            vals.insert(name.clone(), v);
        }
        let (ro, rn) = Interp::new(m, vals.clone()).step()?;
        let (no, nn) = netlist_step(n, m, &vals);
        if ro != no || rn != nn {
            return Ok(Equivalence {
                points: p + 1,
                exhaustive,
                mismatch: Some(format!("at {vals:?}: RTL gives outputs {ro:?}, next {rn:?}; netlist {no:?}, next {nn:?}")),
            });
        }
    }
    Ok(Equivalence { points, exhaustive, mismatch: None })
}
