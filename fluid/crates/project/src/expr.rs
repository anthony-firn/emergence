//! Unit-checked interval expressions.
//!
//! The glue between disciplines: `q1.power * (q1.rjc + tim.r + sink.r)`.
//! Literals carry units in brackets (`0.1[mm]`, `3[W/(m*K)]`), identifiers
//! are node ids (dots allowed), and a call whose name is a node id evaluates a
//! tabulated property at its argument (`mat.air.k(ambient)`). Adding watts to
//! kelvin is an error at evaluation, not a silently wrong number.

use crate::units::{parse_unit, Dim};
use crate::value::Val;

#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Num(Val),
    Var(String),
    Neg(Box<Expr>),
    Bin(char, Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, i32),
    Call(String, Vec<Expr>),
}

/// What an expression needs from its surroundings.
pub trait Scope {
    /// The value of node `id`.
    fn var(&self, id: &str) -> Result<Val, String>;
    /// A tabulated function `id(arg)`, e.g. a material property against
    /// temperature.
    fn table(&self, id: &str, arg: &Val) -> Result<Val, String>;
}

const FUNCTIONS: &[&str] = &["sqrt", "exp", "ln", "tanh", "abs", "min", "max"];

impl Expr {
    pub fn parse(text: &str) -> Result<Expr, String> {
        let toks = tokenize(text)?;
        let mut p = Parser { t: &toks, i: 0 };
        let e = p.sum()?;
        if p.i != toks.len() {
            return Err(format!("`{text}`: unexpected {:?}", toks[p.i]));
        }
        Ok(e)
    }

    /// Node ids this expression reads (variables and tables), in order of
    /// first appearance.
    pub fn deps(&self) -> Vec<String> {
        let mut out = Vec::new();
        self.collect(&mut out);
        out
    }

    fn collect(&self, out: &mut Vec<String>) {
        let mut push = |s: &String| {
            if !out.contains(s) {
                out.push(s.clone());
            }
        };
        match self {
            Expr::Num(_) => {}
            Expr::Var(v) => push(v),
            Expr::Neg(e) | Expr::Pow(e, _) => e.collect(out),
            Expr::Bin(_, a, b) => {
                a.collect(out);
                b.collect(out);
            }
            Expr::Call(f, args) => {
                if !FUNCTIONS.contains(&f.as_str()) {
                    push(f);
                }
                for a in args {
                    a.collect(out);
                }
            }
        }
    }

    pub fn eval(&self, s: &dyn Scope) -> Result<Val, String> {
        match self {
            Expr::Num(v) => Ok(*v),
            Expr::Var(v) => s.var(v),
            Expr::Neg(e) => Ok(e.eval(s)?.neg()),
            Expr::Pow(e, n) => e.eval(s)?.powi(*n),
            Expr::Bin(op, a, b) => {
                let (a, b) = (a.eval(s)?, b.eval(s)?);
                match op {
                    '+' => a.add(&b),
                    '-' => a.sub(&b),
                    '*' => Ok(a.mul(&b)),
                    '/' => a.div(&b),
                    _ => unreachable!(),
                }
            }
            Expr::Call(f, args) => {
                let vals = args.iter().map(|a| a.eval(s)).collect::<Result<Vec<_>, _>>()?;
                let one = || -> Result<&Val, String> {
                    match vals.as_slice() {
                        [v] => Ok(v),
                        _ => Err(format!("{f}() takes one argument")),
                    }
                };
                match f.as_str() {
                    "sqrt" => one()?.sqrt(),
                    "exp" => one()?.exp(),
                    "ln" => one()?.ln(),
                    "tanh" => one()?.tanh(),
                    "abs" => Ok(one()?.abs()),
                    "min" | "max" => {
                        let mut it = vals.iter();
                        let mut acc = *it.next().ok_or_else(|| format!("{f}() needs arguments"))?;
                        for v in it {
                            acc = if f == "min" { acc.min(v)? } else { acc.max(v)? };
                        }
                        Ok(acc)
                    }
                    _ => s.table(f, one()?),
                }
            }
        }
    }

    /// Rewrite every reference to node `from` as `to` (an id rename).
    pub fn rename(&mut self, from: &str, to: &str) {
        match self {
            Expr::Num(_) => {}
            Expr::Var(v) => {
                if v == from {
                    *v = to.into();
                }
            }
            Expr::Neg(e) | Expr::Pow(e, _) => e.rename(from, to),
            Expr::Bin(_, a, b) => {
                a.rename(from, to);
                b.rename(from, to);
            }
            Expr::Call(f, args) => {
                if f == from {
                    *f = to.into();
                }
                for a in args {
                    a.rename(from, to);
                }
            }
        }
    }
}

impl std::fmt::Display for Expr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Expr::Num(v) => {
                if v.dim.is_none() {
                    write!(f, "{}", v.nom)
                } else {
                    write!(f, "{}[{}]", v.nom, v.dim.symbol())
                }
            }
            Expr::Var(v) => write!(f, "{v}"),
            Expr::Neg(e) => write!(f, "-({e})"),
            Expr::Pow(e, n) => write!(f, "({e})^{n}"),
            Expr::Bin(op, a, b) => write!(f, "({a} {op} {b})"),
            Expr::Call(name, args) => {
                write!(f, "{name}(")?;
                for (i, a) in args.iter().enumerate() {
                    if i > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{a}")?;
                }
                write!(f, ")")
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Num(f64),
    Unit(String),
    Id(String),
    Op(char),
}

fn tokenize(text: &str) -> Result<Vec<Tok>, String> {
    let c: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < c.len() {
        let ch = c[i];
        if ch.is_whitespace() {
            i += 1;
        } else if ch.is_ascii_digit() || (ch == '.' && c.get(i + 1).is_some_and(|n| n.is_ascii_digit())) {
            let start = i;
            while i < c.len() && (c[i].is_ascii_digit() || c[i] == '.') {
                i += 1;
            }
            if i < c.len() && (c[i] == 'e' || c[i] == 'E') {
                let save = i;
                i += 1;
                if i < c.len() && (c[i] == '-' || c[i] == '+') {
                    i += 1;
                }
                if i < c.len() && c[i].is_ascii_digit() {
                    while i < c.len() && c[i].is_ascii_digit() {
                        i += 1;
                    }
                } else {
                    i = save;
                }
            }
            let s: String = c[start..i].iter().collect();
            out.push(Tok::Num(s.parse().map_err(|_| format!("bad number `{s}`"))?));
        } else if ch == '[' {
            let start = i + 1;
            while i < c.len() && c[i] != ']' {
                i += 1;
            }
            if i == c.len() {
                return Err(format!("`{text}`: missing `]`"));
            }
            out.push(Tok::Unit(c[start..i].iter().collect()));
            i += 1;
        } else if ch.is_alphabetic() || ch == '_' {
            let start = i;
            while i < c.len() && (c[i].is_alphanumeric() || c[i] == '_' || c[i] == '.') {
                i += 1;
            }
            out.push(Tok::Id(c[start..i].iter().collect()));
        } else if "+-*/^(),".contains(ch) {
            out.push(Tok::Op(ch));
            i += 1;
        } else {
            return Err(format!("`{text}`: unexpected character `{ch}`"));
        }
    }
    Ok(out)
}

struct Parser<'a> {
    t: &'a [Tok],
    i: usize,
}

impl Parser<'_> {
    fn op(&self) -> Option<char> {
        match self.t.get(self.i) {
            Some(Tok::Op(c)) => Some(*c),
            _ => None,
        }
    }
    fn expect(&mut self, c: char) -> Result<(), String> {
        if self.op() == Some(c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("expected `{c}`"))
        }
    }
    fn sum(&mut self) -> Result<Expr, String> {
        let mut e = self.product()?;
        while let Some(op @ ('+' | '-')) = self.op() {
            self.i += 1;
            e = Expr::Bin(op, Box::new(e), Box::new(self.product()?));
        }
        Ok(e)
    }
    fn product(&mut self) -> Result<Expr, String> {
        let mut e = self.unary()?;
        while let Some(op @ ('*' | '/')) = self.op() {
            self.i += 1;
            e = Expr::Bin(op, Box::new(e), Box::new(self.unary()?));
        }
        Ok(e)
    }
    fn unary(&mut self) -> Result<Expr, String> {
        if self.op() == Some('-') {
            self.i += 1;
            return Ok(Expr::Neg(Box::new(self.unary()?)));
        }
        if self.op() == Some('+') {
            self.i += 1;
        }
        self.power()
    }
    fn power(&mut self) -> Result<Expr, String> {
        let base = self.atom()?;
        if self.op() == Some('^') {
            self.i += 1;
            let mut sign = 1;
            let paren = self.op() == Some('(');
            if paren {
                self.i += 1;
            }
            if self.op() == Some('-') {
                sign = -1;
                self.i += 1;
            }
            let n = match self.t.get(self.i) {
                Some(Tok::Num(n)) if n.fract() == 0.0 => *n as i32,
                _ => return Err("exponents must be integers (use sqrt for halves)".into()),
            };
            self.i += 1;
            if paren {
                self.expect(')')?;
            }
            return Ok(Expr::Pow(Box::new(base), sign * n));
        }
        Ok(base)
    }
    fn atom(&mut self) -> Result<Expr, String> {
        match self.t.get(self.i).cloned() {
            Some(Tok::Num(n)) => {
                self.i += 1;
                if let Some(Tok::Unit(u)) = self.t.get(self.i) {
                    self.i += 1;
                    let unit = parse_unit(u)?;
                    if unit.is_affine() {
                        return Err(format!(
                            "`{n}[{u}]`: affine units are not allowed inside expressions; write {}[K]",
                            unit.to_si(n)
                        ));
                    }
                    return Ok(Expr::Num(Val::exact(unit.to_si(n), unit.dim)));
                }
                Ok(Expr::Num(Val::exact(n, Dim::NONE)))
            }
            Some(Tok::Id(id)) => {
                self.i += 1;
                if id == "pi" {
                    return Ok(Expr::Num(Val::exact(std::f64::consts::PI, Dim::NONE)));
                }
                if self.op() == Some('(') {
                    self.i += 1;
                    let mut args = Vec::new();
                    if self.op() != Some(')') {
                        loop {
                            args.push(self.sum()?);
                            if self.op() == Some(',') {
                                self.i += 1;
                            } else {
                                break;
                            }
                        }
                    }
                    self.expect(')')?;
                    return Ok(Expr::Call(id, args));
                }
                Ok(Expr::Var(id))
            }
            Some(Tok::Op('(')) => {
                self.i += 1;
                let e = self.sum()?;
                self.expect(')')?;
                Ok(e)
            }
            Some(Tok::Unit(u)) => Err(format!("unit `[{u}]` must follow a number")),
            other => Err(format!("expected a value, found {other:?}")),
        }
    }
}
