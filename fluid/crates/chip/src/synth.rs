//! Logic synthesis: RTL to a gate netlist of the characterised cells.
//!
//! Every expression is blasted to bits in a hash-consed graph of six node
//! kinds (constant, input, NOT, NAND2, NOR2, XOR2, and the full adder's sum
//! and carry), simplified as it is built: constants fold, double inversions
//! cancel, `x & x = x`, `x ^ x = 0`, identical gates are shared. Additions
//! and comparisons become ripple chains of full-adder cells; multiplexers
//! are three NAND2s and an inverter per bit; `case` is a priority chain of
//! equality tests. Each register bit becomes a D flip-flop whose D input is
//! its next-state function. The result maps one-to-one onto the library:
//! the cells the netlist names are the cells that were simulated at
//! transistor level and characterised.

use std::collections::{BTreeMap, HashMap};

use semiconductor::gates::CellKind;
use semiconductor::timing::{GateNet, Inst};

use crate::rtl::{Dir, Expr, Module, Stmt};

pub type L = u32;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum Node {
    Const(bool),
    In(u32),
    Not(L),
    Nand(L, L),
    Nor(L, L),
    Xor(L, L),
    FaSum(L, L, L),
    FaCarry(L, L, L),
}

/// The bit-level logic graph under construction.
#[derive(Default)]
struct Graph {
    nodes: Vec<Node>,
    index: HashMap<Node, L>,
}

impl Graph {
    fn add(&mut self, n: Node) -> L {
        if let Some(&i) = self.index.get(&n) {
            return i;
        }
        let i = self.nodes.len() as L;
        self.nodes.push(n);
        self.index.insert(n, i);
        i
    }
    fn konst(&mut self, b: bool) -> L {
        self.add(Node::Const(b))
    }
    fn as_const(&self, a: L) -> Option<bool> {
        match self.nodes[a as usize] {
            Node::Const(b) => Some(b),
            _ => None,
        }
    }
    fn is_not_of(&self, a: L, b: L) -> bool {
        self.nodes[a as usize] == Node::Not(b) || self.nodes[b as usize] == Node::Not(a)
    }
    fn not(&mut self, a: L) -> L {
        match self.nodes[a as usize] {
            Node::Const(b) => self.konst(!b),
            Node::Not(x) => x,
            _ => self.add(Node::Not(a)),
        }
    }
    fn nand(&mut self, a: L, b: L) -> L {
        match (self.as_const(a), self.as_const(b)) {
            (Some(false), _) | (_, Some(false)) => return self.konst(true),
            (Some(true), _) => return self.not(b),
            (_, Some(true)) => return self.not(a),
            _ => {}
        }
        if a == b {
            return self.not(a);
        }
        if self.is_not_of(a, b) {
            return self.konst(true);
        }
        let (a, b) = (a.min(b), a.max(b));
        self.add(Node::Nand(a, b))
    }
    fn nor(&mut self, a: L, b: L) -> L {
        match (self.as_const(a), self.as_const(b)) {
            (Some(true), _) | (_, Some(true)) => return self.konst(false),
            (Some(false), _) => return self.not(b),
            (_, Some(false)) => return self.not(a),
            _ => {}
        }
        if a == b {
            return self.not(a);
        }
        if self.is_not_of(a, b) {
            return self.konst(false);
        }
        let (a, b) = (a.min(b), a.max(b));
        self.add(Node::Nor(a, b))
    }
    fn and(&mut self, a: L, b: L) -> L {
        let n = self.nand(a, b);
        self.not(n)
    }
    fn or(&mut self, a: L, b: L) -> L {
        let n = self.nor(a, b);
        self.not(n)
    }
    fn xor(&mut self, a: L, b: L) -> L {
        match (self.as_const(a), self.as_const(b)) {
            (Some(x), Some(y)) => return self.konst(x ^ y),
            (Some(false), _) => return b,
            (_, Some(false)) => return a,
            (Some(true), _) => return self.not(b),
            (_, Some(true)) => return self.not(a),
            _ => {}
        }
        if a == b {
            return self.konst(false);
        }
        if self.is_not_of(a, b) {
            return self.konst(true);
        }
        let (a, b) = (a.min(b), a.max(b));
        self.add(Node::Xor(a, b))
    }
    /// Full adder: (sum, carry). Constant inputs reduce to half adders and
    /// wires.
    fn fa(&mut self, a: L, b: L, c: L) -> (L, L) {
        let mut v = [a, b, c];
        v.sort();
        let consts: Vec<bool> = v.iter().filter_map(|&x| self.as_const(x)).collect();
        let vars: Vec<L> = v.iter().copied().filter(|&x| self.as_const(x).is_none()).collect();
        let ones = consts.iter().filter(|b| **b).count();
        match vars.len() {
            0 => {
                let s = consts.iter().fold(false, |x, y| x ^ y);
                let n = consts.iter().filter(|b| **b).count();
                (self.konst(s), self.konst(n >= 2))
            }
            1 => {
                let x = vars[0];
                match ones {
                    0 => (x, self.konst(false)),
                    1 => (self.not(x), x),
                    _ => (x, self.konst(true)),
                }
            }
            2 => {
                let (x, y) = (vars[0], vars[1]);
                if ones == 0 {
                    (self.xor(x, y), self.and(x, y))
                } else {
                    let s = self.xor(x, y);
                    (self.not(s), self.or(x, y))
                }
            }
            _ => {
                let (a, b, c) = (v[0], v[1], v[2]);
                (self.add(Node::FaSum(a, b, c)), self.add(Node::FaCarry(a, b, c)))
            }
        }
    }
    fn mux(&mut self, s: L, a: L, b: L) -> L {
        // s ? a : b
        match self.as_const(s) {
            Some(true) => return a,
            Some(false) => return b,
            None => {}
        }
        if a == b {
            return a;
        }
        let ns = self.not(s);
        let x = self.nand(s, a);
        let y = self.nand(ns, b);
        self.nand(x, y)
    }
}

/// Cone resynthesis.
///
/// Every node whose function depends on at most `MAX_SUPPORT` inputs (input
/// bits, register outputs) has its truth table computed bit-parallel (one
/// 64-bit word), and is rebuilt from it by Shannon expansion into a balanced
/// tree: equal cofactors drop the variable, complementary ones become an
/// XOR, a constant cofactor becomes an AND or an OR, otherwise a
/// multiplexer; identical subfunctions are shared. The rebuilt cone is used
/// where it is shallower by the cells' characterised delays. A `?:` chain
/// decoding a 4-bit counter, 13 multiplexers deep, becomes 4 levels.
pub const MAX_SUPPORT: usize = 6;

/// Rough per-cell delay (ps) at 1.2 V for choosing between implementations
/// (characterised values, analysis/gates.md: INV 40, NAND2 52, NOR2 50,
/// XOR2 87, FA 177).
fn cell_delay(n: &Node) -> f64 {
    match n {
        Node::Not(_) => 40.0,
        Node::Nand(..) => 52.0,
        Node::Nor(..) => 50.0,
        Node::Xor(..) => 87.0,
        Node::FaSum(..) | Node::FaCarry(..) => 177.0,
        _ => 0.0,
    }
}

fn children(n: &Node) -> Vec<L> {
    match *n {
        Node::Not(a) => vec![a],
        Node::Nand(a, b) | Node::Nor(a, b) | Node::Xor(a, b) => vec![a, b],
        Node::FaSum(a, b, c) | Node::FaCarry(a, b, c) => vec![a, b, c],
        _ => vec![],
    }
}

fn full(k: usize) -> u64 {
    if k >= 6 { u64::MAX } else { (1u64 << (1usize << k)) - 1 }
}

/// Masks of the variables over 2^k rows.
const VAR: [u64; 6] = [
    0xAAAA_AAAA_AAAA_AAAA,
    0xCCCC_CCCC_CCCC_CCCC,
    0xF0F0_F0F0_F0F0_F0F0,
    0xFF00_FF00_FF00_FF00,
    0xFFFF_0000_FFFF_0000,
    0xFFFF_FFFF_0000_0000,
];

/// Cofactors of `tt` on variable `v` (of `k`), as tables over k−1 vars
/// with `v` removed (the higher variables shift down).
fn cofactors(tt: u64, v: usize, k: usize) -> (u64, u64) {
    let rows = 1usize << k;
    let (mut f0, mut f1) = (0u64, 0u64);
    let mut j = 0;
    for r in 0..rows {
        if (r >> v) & 1 == 0 {
            let r1 = r | (1 << v);
            f0 |= ((tt >> r) & 1) << j;
            f1 |= ((tt >> r1) & 1) << j;
            j += 1;
        }
    }
    (f0, f1)
}

impl Graph {
    fn depth(&self, memo: &mut Vec<f64>) {
        memo.clear();
        for n in &self.nodes {
            let d = children(n).iter().map(|&c| memo[c as usize]).fold(0.0, f64::max) + cell_delay(n);
            memo.push(d);
        }
    }

    /// Build the function `tt` of `vars` (graph nodes), sharing subtables.
    fn from_table(&mut self, tt: u64, vars: &[L], memo: &mut HashMap<(u64, usize), L>) -> L {
        let k = vars.len();
        let tt = tt & full(k);
        if tt == 0 {
            return self.konst(false);
        }
        if tt == full(k) {
            return self.konst(true);
        }
        if let Some(&n) = memo.get(&(tt, k)) {
            return n;
        }
        let v = k - 1;
        let x = vars[v];
        let (f0, f1) = cofactors(tt, v, k);
        let rest = &vars[..v];
        let m1 = full(v);
        let r = if f0 == f1 {
            self.from_table(f0, rest, memo)
        } else if f0 == !f1 & m1 {
            let a = self.from_table(f0, rest, memo);
            self.xor(x, a)
        } else if f0 == 0 {
            let a = self.from_table(f1, rest, memo);
            self.and(x, a)
        } else if f1 == 0 {
            let a = self.from_table(f0, rest, memo);
            let nx = self.not(x);
            self.and(nx, a)
        } else if f0 == m1 {
            let a = self.from_table(f1, rest, memo);
            let nx = self.not(x);
            self.or(nx, a)
        } else if f1 == m1 {
            let a = self.from_table(f0, rest, memo);
            self.or(x, a)
        } else {
            let a = self.from_table(f1, rest, memo);
            let b = self.from_table(f0, rest, memo);
            self.mux(x, a, b)
        };
        memo.insert((tt, k), r);
        r
    }
}

/// Rebuild `g` with small-support cones resynthesized where that is
/// shallower. Returns the new graph and the map from old to new nodes.
fn resynthesize(g: &Graph) -> (Graph, Vec<L>) {
    // Supports (sorted input ids), capped one past the limit.
    let mut support: Vec<Vec<u32>> = Vec::with_capacity(g.nodes.len());
    for n in &g.nodes {
        let s = match *n {
            Node::In(k) => vec![k],
            Node::Const(_) => vec![],
            _ => {
                let mut s: Vec<u32> = children(n).iter().flat_map(|&c| support[c as usize].iter().copied()).collect();
                s.sort();
                s.dedup();
                s.truncate(MAX_SUPPORT + 1);
                s
            }
        };
        support.push(s);
    }
    let mut depth = Vec::new();
    g.depth(&mut depth);
    let mut out = Graph::default();
    let mut map: Vec<L> = Vec::with_capacity(g.nodes.len());
    let mut tts: HashMap<L, u64> = HashMap::new();
    let mut in_node: HashMap<u32, L> = HashMap::new();
    for (i, n) in g.nodes.iter().enumerate() {
        let built = match *n {
            Node::In(k) => {
                let x = out.add(Node::In(k));
                in_node.insert(k, x);
                x
            }
            Node::Const(b) => out.konst(b),
            Node::Not(a) => out.not(map[a as usize]),
            Node::Nand(a, b) => out.nand(map[a as usize], map[b as usize]),
            Node::Nor(a, b) => out.nor(map[a as usize], map[b as usize]),
            Node::Xor(a, b) => out.xor(map[a as usize], map[b as usize]),
            Node::FaSum(a, b, c) => out.fa(map[a as usize], map[b as usize], map[c as usize]).0,
            Node::FaCarry(a, b, c) => out.fa(map[a as usize], map[b as usize], map[c as usize]).1,
        };
        let sup = &support[i];
        let mut chosen = built;
        if sup.len() <= MAX_SUPPORT && !matches!(n, Node::In(_) | Node::Const(_)) {
            // Truth table of the original node over its support.
            let pos = |k: u32| sup.iter().position(|&s| s == k).unwrap();
            let tt = {
                let mut memo: HashMap<L, u64> = HashMap::new();
                fn eval(g: &Graph, x: L, pos: &dyn Fn(u32) -> usize, memo: &mut HashMap<L, u64>) -> u64 {
                    if let Some(&t) = memo.get(&x) {
                        return t;
                    }
                    let t = match g.nodes[x as usize] {
                        Node::In(k) => VAR[pos(k)],
                        Node::Const(b) => if b { u64::MAX } else { 0 },
                        Node::Not(a) => !eval(g, a, pos, memo),
                        Node::Nand(a, b) => !(eval(g, a, pos, memo) & eval(g, b, pos, memo)),
                        Node::Nor(a, b) => !(eval(g, a, pos, memo) | eval(g, b, pos, memo)),
                        Node::Xor(a, b) => eval(g, a, pos, memo) ^ eval(g, b, pos, memo),
                        Node::FaSum(a, b, c) => eval(g, a, pos, memo) ^ eval(g, b, pos, memo) ^ eval(g, c, pos, memo),
                        Node::FaCarry(a, b, c) => {
                            let (x, y, z) = (eval(g, a, pos, memo), eval(g, b, pos, memo), eval(g, c, pos, memo));
                            (x & y) | (x & z) | (y & z)
                        }
                    };
                    memo.insert(x, t);
                    t
                }
                eval(g, i as L, &pos, &mut memo) & full(sup.len())
            };
            let vars: Vec<L> = sup.iter().map(|k| in_node[k]).collect();
            let mut memo = HashMap::new();
            let alt = out.from_table(tt, &vars, &mut memo);
            if alt != built {
                let mut d = Vec::new();
                out.depth(&mut d);
                if d[alt as usize] + 1e-9 < d[built as usize] {
                    chosen = alt;
                }
            }
            tts.insert(i as L, tt);
        }
        map.push(chosen);
    }
    let _ = depth;
    (out, map)
}

type Bits = Vec<L>;

/// Width of an expression by Verilog's self-determined rules (unsized
/// constants take the bits they need).
pub fn width(m: &Module, e: &Expr) -> Result<usize, String> {
    Ok(match e {
        Expr::Const(_, Some(w)) => *w,
        Expr::Const(v, None) => (128 - v.leading_zeros() as usize).max(1),
        Expr::Ident(n) => *m.widths.get(n).ok_or(format!("unknown signal `{n}`"))?,
        Expr::Slice(_, a, b) => a - b + 1,
        Expr::Concat(p) => p.iter().map(|x| width(m, x)).sum::<Result<usize, String>>()?,
        Expr::Repeat(n, x) => n * width(m, x)?,
        Expr::Unary(op, a) => match *op {
            "~" | "-" => width(m, a)?,
            _ => 1,
        },
        Expr::Binary(op, a, b) => match *op {
            "&" | "|" | "^" | "+" | "-" | "*" => width(m, a)?.max(width(m, b)?),
            "<<" | ">>" => width(m, a)?,
            _ => 1,
        },
        Expr::Cond(_, a, b) => width(m, a)?.max(width(m, b)?),
    })
}

struct Ctx<'a> {
    m: &'a Module,
    g: Graph,
    /// Bits of inputs and register outputs.
    sources: BTreeMap<String, Bits>,
    /// Bit drivers of wires: (assign index, bit within its value) per bit.
    drivers: BTreeMap<String, Vec<Option<(usize, usize)>>>,
    memo: HashMap<usize, Bits>,
    visiting: Vec<usize>,
}

fn fit(g: &mut Graph, mut v: Bits, w: usize) -> Bits {
    let z = g.konst(false);
    v.resize(w, z);
    v
}

impl Ctx<'_> {
    fn signal(&mut self, n: &str) -> Result<Bits, String> {
        if let Some(b) = self.sources.get(n) {
            return Ok(b.clone());
        }
        let drv = self.drivers.get(n).cloned().ok_or(format!("`{n}` is used but never assigned"))?;
        let mut out = Vec::with_capacity(drv.len());
        for (bit, d) in drv.iter().enumerate() {
            let (ai, off) = d.ok_or(format!("bit {bit} of `{n}` is never assigned"))?;
            let v = self.assign_value(ai)?;
            out.push(v[off]);
        }
        Ok(out)
    }

    fn assign_value(&mut self, ai: usize) -> Result<Bits, String> {
        if let Some(v) = self.memo.get(&ai) {
            return Ok(v.clone());
        }
        if self.visiting.contains(&ai) {
            return Err(format!("combinational loop through `{}`", self.m.assigns[ai].0));
        }
        self.visiting.push(ai);
        let (name, part, e) = self.m.assigns[ai].clone();
        let w = match part {
            Some((a, b)) => a - b + 1,
            None => self.m.widths[&name],
        };
        let ctx = w.max(width(self.m, &e)?);
        let v = self.bits(&e, ctx)?;
        let v = v[..w].to_vec();
        self.visiting.pop();
        self.memo.insert(ai, v.clone());
        Ok(v)
    }

    /// `e` evaluated at context width `w` (exactly `w` bits, LSB first).
    fn bits(&mut self, e: &Expr, w: usize) -> Result<Bits, String> {
        let m = self.m;
        let v: Bits = match e {
            Expr::Const(c, _) => (0..w).map(|i| self.g.konst(i < 128 && (c >> i) & 1 == 1)).collect(),
            Expr::Ident(n) => {
                let s = self.signal(n)?;
                fit(&mut self.g, s, w)
            }
            Expr::Slice(x, msb, lsb) => {
                let sw = width(m, x)?;
                if *msb >= sw {
                    return Err(format!("select [{msb}:{lsb}] is outside a {sw}-bit value"));
                }
                let s = self.bits(x, sw)?;
                fit(&mut self.g, s[*lsb..=*msb].to_vec(), w)
            }
            Expr::Concat(parts) => {
                let mut out = Vec::new();
                for p in parts.iter().rev() {
                    let pw = width(m, p)?;
                    out.extend(self.bits(p, pw)?);
                }
                fit(&mut self.g, out, w)
            }
            Expr::Repeat(n, x) => {
                let xw = width(m, x)?;
                let one = self.bits(x, xw)?;
                let mut out = Vec::new();
                for _ in 0..*n {
                    out.extend(one.iter().copied());
                }
                fit(&mut self.g, out, w)
            }
            Expr::Unary(op, a) => match *op {
                "~" => {
                    let a = self.bits(a, w)?;
                    a.into_iter().map(|x| self.g.not(x)).collect()
                }
                "-" => {
                    let a = self.bits(a, w)?;
                    let na: Bits = a.into_iter().map(|x| self.g.not(x)).collect();
                    let zero = vec![self.g.konst(false); w];
                    let one = self.g.konst(true);
                    self.add(&na, &zero, one).0
                }
                "!" | "|" => {
                    let aw = width(m, a)?;
                    let a = self.bits(a, aw)?;
                    let r = self.reduce_or(&a);
                    let r = if *op == "!" { self.g.not(r) } else { r };
                    fit(&mut self.g, vec![r], w)
                }
                "&" => {
                    let aw = width(m, a)?;
                    let a = self.bits(a, aw)?;
                    let r = a.iter().skip(1).fold(a[0], |acc, &x| self.g.and(acc, x));
                    fit(&mut self.g, vec![r], w)
                }
                "^" => {
                    let aw = width(m, a)?;
                    let a = self.bits(a, aw)?;
                    let r = a.iter().skip(1).fold(a[0], |acc, &x| self.g.xor(acc, x));
                    fit(&mut self.g, vec![r], w)
                }
                _ => unreachable!(),
            },
            Expr::Binary(op, a, b) => match *op {
                "&" | "|" | "^" => {
                    let (a, b) = (self.bits(a, w)?, self.bits(b, w)?);
                    a.into_iter()
                        .zip(b)
                        .map(|(x, y)| match *op {
                            "&" => self.g.and(x, y),
                            "|" => self.g.or(x, y),
                            _ => self.g.xor(x, y),
                        })
                        .collect()
                }
                "+" => {
                    let (a, b) = (self.bits(a, w)?, self.bits(b, w)?);
                    let z = self.g.konst(false);
                    self.add(&a, &b, z).0
                }
                "-" => {
                    let (a, b) = (self.bits(a, w)?, self.bits(b, w)?);
                    let nb: Bits = b.into_iter().map(|x| self.g.not(x)).collect();
                    let one = self.g.konst(true);
                    self.add(&a, &nb, one).0
                }
                "*" => {
                    let (a, b) = (self.bits(a, w)?, self.bits(b, w)?);
                    let mut acc = vec![self.g.konst(false); w];
                    for (i, &bi) in b.iter().enumerate() {
                        if self.g.as_const(bi) == Some(false) {
                            continue;
                        }
                        let z = self.g.konst(false);
                        let mut pp = vec![z; i];
                        pp.extend(a.iter().take(w - i).map(|&x| self.g.and(x, bi)));
                        let z = self.g.konst(false);
                        acc = self.add(&acc, &pp, z).0;
                    }
                    acc
                }
                "==" | "!=" => {
                    let ow = width(m, a)?.max(width(m, b)?);
                    let (a, b) = (self.bits(a, ow)?, self.bits(b, ow)?);
                    let diff: Bits = a.into_iter().zip(b).map(|(x, y)| self.g.xor(x, y)).collect();
                    let any = self.reduce_or(&diff);
                    let r = if *op == "==" { self.g.not(any) } else { any };
                    fit(&mut self.g, vec![r], w)
                }
                "<" | "<=" | ">" | ">=" => {
                    let ow = width(m, a)?.max(width(m, b)?);
                    let (a, b) = (self.bits(a, ow)?, self.bits(b, ow)?);
                    // a >= b  ⇔  carry out of a + ~b + 1
                    let (x, y) = if *op == "<" || *op == ">=" { (a, b) } else { (b, a) };
                    let ny: Bits = y.into_iter().map(|v| self.g.not(v)).collect();
                    let one = self.g.konst(true);
                    let (_, ge) = self.add(&x, &ny, one);
                    let r = if *op == "<" || *op == ">" { self.g.not(ge) } else { ge };
                    fit(&mut self.g, vec![r], w)
                }
                "&&" | "||" => {
                    let (aw, bw) = (width(m, a)?, width(m, b)?);
                    let (a, b) = (self.bits(a, aw)?, self.bits(b, bw)?);
                    let (x, y) = (self.reduce_or(&a), self.reduce_or(&b));
                    let r = if *op == "&&" { self.g.and(x, y) } else { self.g.or(x, y) };
                    fit(&mut self.g, vec![r], w)
                }
                "<<" | ">>" => {
                    let a = self.bits(a, w)?;
                    let bw = width(m, b)?;
                    let sh = self.bits(b, bw)?;
                    // Barrel shifter, one stage per shift-amount bit.
                    let mut cur = a;
                    for (k, &s) in sh.iter().enumerate() {
                        let d = 1usize << k.min(63);
                        let z = self.g.konst(false);
                        let shifted: Bits = (0..w)
                            .map(|i| {
                                if *op == "<<" {
                                    if i >= d { cur[i - d] } else { z }
                                } else if i + d < w {
                                    cur[i + d]
                                } else {
                                    z
                                }
                            })
                            .collect();
                        cur = cur.iter().zip(&shifted).map(|(&x, &y)| self.g.mux(s, y, x)).collect();
                    }
                    cur
                }
                other => return Err(format!("operator `{other}` is not supported")),
            },
            Expr::Cond(c, a, b) => {
                let cw = width(m, c)?;
                let c = self.bits(c, cw)?;
                let s = self.reduce_or(&c);
                let (a, b) = (self.bits(a, w)?, self.bits(b, w)?);
                a.into_iter().zip(b).map(|(x, y)| self.g.mux(s, x, y)).collect()
            }
        };
        debug_assert_eq!(v.len(), w);
        Ok(v)
    }

    fn reduce_or(&mut self, a: &[L]) -> L {
        let mut acc = a[0];
        for &x in &a[1..] {
            acc = self.g.or(acc, x);
        }
        acc
    }

    fn add(&mut self, a: &[L], b: &[L], cin: L) -> (Bits, L) {
        let mut c = cin;
        let mut s = Vec::with_capacity(a.len());
        for (&x, &y) in a.iter().zip(b) {
            let (sum, carry) = self.g.fa(x, y, c);
            s.push(sum);
            c = carry;
        }
        (s, c)
    }

    /// Next state of the registers through a list of statements.
    fn exec(&mut self, stmts: &[Stmt], next: &mut BTreeMap<String, Bits>) -> Result<(), String> {
        for st in stmts {
            match st {
                Stmt::Assign(n, part, e) => {
                    let total = self.m.widths[n];
                    let (msb, lsb) = part.unwrap_or((total - 1, 0));
                    if msb >= total {
                        return Err(format!("`{n}[{msb}]` is outside its {total} bits"));
                    }
                    let w = msb - lsb + 1;
                    let ctx = w.max(width(self.m, e)?);
                    let v = self.bits(e, ctx)?;
                    let slot = next.get_mut(n).unwrap();
                    slot[lsb..=msb].copy_from_slice(&v[..w]);
                }
                Stmt::If(c, a, b) => {
                    let cw = width(self.m, c)?;
                    let cb = self.bits(c, cw)?;
                    let s = self.reduce_or(&cb);
                    let mut na = next.clone();
                    self.exec(a, &mut na)?;
                    let mut nb = next.clone();
                    self.exec(b, &mut nb)?;
                    for (k, v) in next.iter_mut() {
                        *v = na[k].iter().zip(&nb[k]).map(|(&x, &y)| self.g.mux(s, x, y)).collect();
                    }
                }
                Stmt::Case(sel, arms, default) => {
                    // A priority chain: the first matching arm wins.
                    let mut chain: Vec<(Expr, Vec<Stmt>)> = Vec::new();
                    for (labels, body) in arms {
                        let mut cond: Option<Expr> = None;
                        for l in labels {
                            let eq = Expr::Binary("==", Box::new(sel.clone()), Box::new(l.clone()));
                            cond = Some(match cond {
                                None => eq,
                                Some(c) => Expr::Binary("||", Box::new(c), Box::new(eq)),
                            });
                        }
                        chain.push((cond.unwrap(), body.clone()));
                    }
                    let mut tail = default.clone();
                    for (c, body) in chain.into_iter().rev() {
                        tail = vec![Stmt::If(c, body, tail)];
                    }
                    self.exec(&tail, next)?;
                }
            }
        }
        Ok(())
    }
}

/// A register bit.
#[derive(Clone, Debug)]
pub struct Flop {
    pub name: String,
    pub d: usize,
    pub q: usize,
}

/// A synthesized design.
#[derive(Clone, Debug)]
pub struct Netlist {
    pub name: String,
    /// Combinational logic. Its primary inputs are, in order: the module's
    /// input bits, the flip-flop outputs, then the constant nets.
    pub comb: GateNet,
    pub inputs: Vec<(String, Vec<usize>)>,
    pub outputs: Vec<(String, Vec<usize>)>,
    pub flops: Vec<Flop>,
    /// Tie-low and tie-high nets (no transitions).
    pub zero: usize,
    pub one: usize,
    pub clock: Option<String>,
}

impl Netlist {
    pub fn cells(&self) -> BTreeMap<&'static str, usize> {
        let mut c = BTreeMap::new();
        for i in &self.comb.insts {
            *c.entry(i.kind.name()).or_insert(0) += 1;
        }
        if !self.flops.is_empty() {
            c.insert("DFF", self.flops.len());
        }
        c
    }
    /// Transistor count: INV 2, NAND2/NOR2 4, XOR2 16, FA 44, DFF 18.
    pub fn transistors(&self) -> usize {
        let comb: usize = self
            .comb
            .insts
            .iter()
            .map(|i| match i.kind {
                CellKind::Inv => 2,
                CellKind::Nand2 | CellKind::Nor2 => 4,
                CellKind::Xor2 => 16,
                CellKind::FullAdder => 44,
            })
            .sum();
        comb + 18 * self.flops.len()
    }
}

/// Synthesize a parsed module.
pub fn synthesize(m: &Module) -> Result<Netlist, String> {
    let mut ctx = Ctx { m, g: Graph::default(), sources: BTreeMap::new(), drivers: BTreeMap::new(), memo: HashMap::new(), visiting: vec![] };
    let mut pi = 0u32;
    let mut input_bits = Vec::new();
    for p in m.ports.iter().filter(|p| p.dir == Dir::Input) {
        if Some(&p.name) == m.clock.as_ref() {
            continue;
        }
        let bits: Bits = (0..p.width).map(|_| {
            let b = ctx.g.add(Node::In(pi));
            pi += 1;
            b
        }).collect();
        input_bits.push((p.name.clone(), bits.clone()));
        ctx.sources.insert(p.name.clone(), bits);
    }
    let mut reg_bits = Vec::new();
    for r in &m.regs {
        let w = m.widths[r];
        let bits: Bits = (0..w).map(|_| {
            let b = ctx.g.add(Node::In(pi));
            pi += 1;
            b
        }).collect();
        reg_bits.push((r.clone(), bits.clone()));
        ctx.sources.insert(r.clone(), bits);
    }
    if let Some(c) = &m.clock {
        if !m.ports.iter().any(|p| &p.name == c && p.dir == Dir::Input && p.width == 1) {
            return Err(format!("clock `{c}` must be a 1-bit input"));
        }
    }
    for (ai, (n, part, _)) in m.assigns.iter().enumerate() {
        let w = m.widths[n];
        let slot = ctx.drivers.entry(n.clone()).or_insert_with(|| vec![None; w]);
        let (msb, lsb) = part.unwrap_or((w - 1, 0));
        if msb >= w {
            return Err(format!("`{n}[{msb}]` is outside its {w} bits"));
        }
        for (k, bit) in (lsb..=msb).enumerate() {
            if slot[bit].is_some() {
                return Err(format!("bit {bit} of `{n}` has two drivers"));
            }
            slot[bit] = Some((ai, k));
        }
    }
    let mut next: BTreeMap<String, Bits> = reg_bits.iter().cloned().collect();
    let clocked = m.clocked.clone();
    ctx.exec(&clocked, &mut next)?;
    let mut out_bits = Vec::new();
    for p in m.ports.iter().filter(|p| p.dir == Dir::Output) {
        out_bits.push((p.name.clone(), ctx.signal(&p.name)?));
    }
    // Resynthesize small cones where that is shallower.
    let (g2, remap) = resynthesize(&ctx.g);
    for (_, b) in out_bits.iter_mut() {
        for x in b.iter_mut() {
            *x = remap[*x as usize];
        }
    }
    for v in next.values_mut() {
        for x in v.iter_mut() {
            *x = remap[*x as usize];
        }
    }
    let reg_bits: Vec<(String, Bits)> = reg_bits.iter().map(|(n, b)| (n.clone(), b.iter().map(|&x| remap[x as usize]).collect())).collect();
    let input_bits: Vec<(String, Bits)> = input_bits.iter().map(|(n, b)| (n.clone(), b.iter().map(|&x| remap[x as usize]).collect())).collect();
    // Emit cells for the logic reachable from the outputs and register
    // inputs, in node order (which is topological).
    let g = &g2;
    let mut used = vec![false; g.nodes.len()];
    let mut stack: Vec<L> = out_bits.iter().flat_map(|(_, b)| b.iter().copied()).collect();
    stack.extend(next.values().flatten().copied());
    while let Some(x) = stack.pop() {
        if used[x as usize] {
            continue;
        }
        used[x as usize] = true;
        match g.nodes[x as usize] {
            Node::Not(a) => stack.push(a),
            Node::Nand(a, b) | Node::Nor(a, b) | Node::Xor(a, b) => stack.extend([a, b]),
            Node::FaSum(a, b, c) | Node::FaCarry(a, b, c) => stack.extend([a, b, c]),
            _ => {}
        }
    }
    let n_pi = pi as usize;
    let zero = n_pi;
    let one = n_pi + 1;
    let mut net_of = vec![usize::MAX; g.nodes.len()];
    let mut n_nets = n_pi + 2;
    let mut insts = Vec::new();
    let mut fa_done: HashMap<(L, L, L), (usize, usize)> = HashMap::new();
    for (i, node) in g.nodes.iter().enumerate() {
        match *node {
            Node::In(k) => net_of[i] = k as usize,
            Node::Const(b) => net_of[i] = if b { one } else { zero },
            _ if !used[i] => {}
            Node::Not(a) => {
                net_of[i] = n_nets;
                insts.push(Inst { kind: CellKind::Inv, ins: vec![net_of[a as usize]], outs: vec![n_nets] });
                n_nets += 1;
            }
            Node::Nand(a, b) | Node::Nor(a, b) | Node::Xor(a, b) => {
                let kind = match node {
                    Node::Nand(..) => CellKind::Nand2,
                    Node::Nor(..) => CellKind::Nor2,
                    _ => CellKind::Xor2,
                };
                net_of[i] = n_nets;
                insts.push(Inst { kind, ins: vec![net_of[a as usize], net_of[b as usize]], outs: vec![n_nets] });
                n_nets += 1;
            }
            Node::FaSum(a, b, c) | Node::FaCarry(a, b, c) => {
                let (s, co) = *fa_done.entry((a, b, c)).or_insert_with(|| {
                    let (s, co) = (n_nets, n_nets + 1);
                    n_nets += 2;
                    insts.push(Inst {
                        kind: CellKind::FullAdder,
                        ins: vec![net_of[a as usize], net_of[b as usize], net_of[c as usize]],
                        outs: vec![s, co],
                    });
                    (s, co)
                });
                net_of[i] = if matches!(node, Node::FaSum(..)) { s } else { co };
            }
        }
    }
    let flops = reg_bits
        .iter()
        .flat_map(|(name, q)| {
            let n = &next[name];
            let w = q.len();
            q.iter().zip(n).enumerate().map(move |(k, (&qb, &db))| (if w == 1 { name.clone() } else { format!("{name}[{k}]") }, qb, db)).collect::<Vec<_>>()
        })
        .map(|(name, qb, db)| Flop { name, q: net_of[qb as usize], d: net_of[db as usize] })
        .collect();
    let map = |bits: &Bits| bits.iter().map(|&b| net_of[b as usize]).collect::<Vec<_>>();
    let mut inputs: Vec<usize> = (0..n_pi).collect();
    inputs.extend([zero, one]);
    Ok(Netlist {
        name: m.name.clone(),
        comb: GateNet { n_nets, inputs, insts },
        inputs: input_bits.iter().map(|(n, b)| (n.clone(), map(b))).collect(),
        outputs: out_bits.iter().map(|(n, b)| (n.clone(), map(b))).collect(),
        flops,
        zero,
        one,
        clock: m.clock.clone(),
    })
}
