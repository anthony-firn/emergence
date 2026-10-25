//! Named parameters and the expression language that makes a model
//! parametric.
//!
//! A [`Params`] table holds named lengths, counts and angles; every number
//! in a JSON model may be written as an expression over them
//! (`"wall + clearance"`, `"pitch * teeth / pi"`). Parameters may be
//! defined in terms of earlier parameters, which is all a design table
//! needs: change one value, rebuild, and every dependent dimension moves.
//!
//! The grammar is the usual precedence-climbing one (Pratt 1973):
//!
//! ```text
//!   expr    := term (('+' | '-') term)*
//!   term    := power (('*' | '/' | '%') power)*
//!   power   := unary ('^' power)?            -- right associative
//!   unary   := ('-' | '+')? primary
//!   primary := number | name | name '(' expr (',' expr)* ')' | '(' expr ')'
//! ```
//!
//! Angles are radians throughout the crate; `rad(x)` converts from degrees
//! and `deg(x)` back.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// A parameter table.
#[derive(Clone, Debug, Default, Deserialize, Serialize)]
#[serde(transparent)]
pub struct Params {
    pub values: BTreeMap<String, f64>,
}

impl Params {
    pub fn new() -> Self {
        Self::default()
    }

    /// Build from `(name, expression)` pairs evaluated in order, so a later
    /// definition may use an earlier one.
    pub fn from_defs(defs: &[(&str, &str)]) -> Result<Self, String> {
        let mut p = Params::new();
        for (name, expr) in defs {
            let v = p.eval(expr)?;
            p.values.insert((*name).to_string(), v);
        }
        Ok(p)
    }

    pub fn set(&mut self, name: &str, value: f64) -> &mut Self {
        self.values.insert(name.to_string(), value);
        self
    }

    /// Set from an expression over the parameters already present.
    pub fn define(&mut self, name: &str, expr: &str) -> Result<&mut Self, String> {
        let v = self.eval(expr)?;
        self.values.insert(name.to_string(), v);
        Ok(self)
    }

    pub fn get(&self, name: &str) -> Option<f64> {
        self.values.get(name).copied()
    }

    /// Evaluate an expression over this table.
    pub fn eval(&self, text: &str) -> Result<f64, String> {
        let tokens = lex(text)?;
        let mut p = Parser { tokens: &tokens, at: 0, params: self };
        let v = p.expr(0)?;
        if p.at != p.tokens.len() {
            return Err(format!("unexpected trailing input in {text:?}"));
        }
        Ok(v)
    }

    /// Evaluate, reporting the expression in the error message.
    pub fn need(&self, text: &str) -> f64 {
        self.eval(text).unwrap_or_else(|e| panic!("{e}"))
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Num(f64),
    Name(String),
    Op(char),
    Open,
    Close,
    Comma,
}

fn lex(text: &str) -> Result<Vec<Tok>, String> {
    let b: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let c = b[i];
        if c.is_whitespace() {
            i += 1;
        } else if c.is_ascii_digit() || c == '.' {
            let start = i;
            while i < b.len() && (b[i].is_ascii_digit() || b[i] == '.') {
                i += 1;
            }
            // Exponent.
            if i < b.len() && (b[i] == 'e' || b[i] == 'E') {
                let save = i;
                i += 1;
                if i < b.len() && (b[i] == '+' || b[i] == '-') {
                    i += 1;
                }
                if i < b.len() && b[i].is_ascii_digit() {
                    while i < b.len() && b[i].is_ascii_digit() {
                        i += 1;
                    }
                } else {
                    i = save;
                }
            }
            let s: String = b[start..i].iter().collect();
            out.push(Tok::Num(s.parse().map_err(|_| format!("bad number {s:?}"))?));
        } else if c.is_alphabetic() || c == '_' {
            let start = i;
            while i < b.len() && (b[i].is_alphanumeric() || b[i] == '_') {
                i += 1;
            }
            out.push(Tok::Name(b[start..i].iter().collect()));
        } else {
            i += 1;
            match c {
                '(' => out.push(Tok::Open),
                ')' => out.push(Tok::Close),
                ',' => out.push(Tok::Comma),
                '+' | '-' | '*' | '/' | '^' | '%' => out.push(Tok::Op(c)),
                _ => return Err(format!("unexpected character {c:?}")),
            }
        }
    }
    Ok(out)
}

struct Parser<'a> {
    tokens: &'a [Tok],
    at: usize,
    params: &'a Params,
}

impl Parser<'_> {
    fn peek(&self) -> Option<&Tok> {
        self.tokens.get(self.at)
    }

    fn expr(&mut self, min_prec: u8) -> Result<f64, String> {
        let mut lhs = self.unary()?;
        while let Some(Tok::Op(op)) = self.peek().cloned() {
            let (prec, right) = match op {
                '+' | '-' => (1, false),
                '*' | '/' | '%' => (2, false),
                '^' => (3, true),
                _ => break,
            };
            if prec < min_prec {
                break;
            }
            self.at += 1;
            let rhs = self.expr(if right { prec } else { prec + 1 })?;
            lhs = match op {
                '+' => lhs + rhs,
                '-' => lhs - rhs,
                '*' => lhs * rhs,
                '/' => lhs / rhs,
                '%' => lhs % rhs,
                '^' => lhs.powf(rhs),
                _ => unreachable!(),
            };
        }
        Ok(lhs)
    }

    fn unary(&mut self) -> Result<f64, String> {
        match self.peek().cloned() {
            Some(Tok::Op('-')) => {
                self.at += 1;
                // Looser than '^', as in mathematical notation: -2^2 = -4.
                Ok(-self.expr(3)?)
            }
            Some(Tok::Op('+')) => {
                self.at += 1;
                self.expr(3)
            }
            _ => self.primary(),
        }
    }

    fn primary(&mut self) -> Result<f64, String> {
        match self.peek().cloned() {
            Some(Tok::Num(v)) => {
                self.at += 1;
                Ok(v)
            }
            Some(Tok::Open) => {
                self.at += 1;
                let v = self.expr(0)?;
                if self.peek() != Some(&Tok::Close) {
                    return Err("missing ')'".into());
                }
                self.at += 1;
                Ok(v)
            }
            Some(Tok::Name(name)) => {
                self.at += 1;
                if self.peek() == Some(&Tok::Open) {
                    self.at += 1;
                    let mut args = Vec::new();
                    if self.peek() != Some(&Tok::Close) {
                        loop {
                            args.push(self.expr(0)?);
                            match self.peek() {
                                Some(Tok::Comma) => self.at += 1,
                                _ => break,
                            }
                        }
                    }
                    if self.peek() != Some(&Tok::Close) {
                        return Err(format!("missing ')' after {name}("));
                    }
                    self.at += 1;
                    call(&name, &args)
                } else {
                    constant(&name)
                        .or_else(|| self.params.get(&name))
                        .ok_or_else(|| format!("unknown parameter {name:?}"))
                }
            }
            other => Err(format!("unexpected {other:?}")),
        }
    }
}

fn constant(name: &str) -> Option<f64> {
    Some(match name {
        "pi" => std::f64::consts::PI,
        "tau" => std::f64::consts::TAU,
        "e" => std::f64::consts::E,
        "inf" => f64::INFINITY,
        _ => return None,
    })
}

fn call(name: &str, a: &[f64]) -> Result<f64, String> {
    let one = |f: fn(f64) -> f64| {
        if a.len() == 1 {
            Ok(f(a[0]))
        } else {
            Err(format!("{name} takes 1 argument"))
        }
    };
    match name {
        "sqrt" => one(f64::sqrt),
        "abs" => one(f64::abs),
        "sin" => one(f64::sin),
        "cos" => one(f64::cos),
        "tan" => one(f64::tan),
        "asin" => one(f64::asin),
        "acos" => one(f64::acos),
        "atan" => one(f64::atan),
        "exp" => one(f64::exp),
        "ln" => one(f64::ln),
        "log10" => one(f64::log10),
        "floor" => one(f64::floor),
        "ceil" => one(f64::ceil),
        "round" => one(f64::round),
        "sign" => one(f64::signum),
        "rad" => one(f64::to_radians),
        "deg" => one(f64::to_degrees),
        "atan2" | "hypot" | "pow" => {
            if a.len() != 2 {
                return Err(format!("{name} takes 2 arguments"));
            }
            Ok(match name {
                "atan2" => a[0].atan2(a[1]),
                "hypot" => a[0].hypot(a[1]),
                _ => a[0].powf(a[1]),
            })
        }
        "min" | "max" => {
            if a.is_empty() {
                return Err(format!("{name} needs arguments"));
            }
            Ok(a.iter().copied().fold(
                if name == "min" { f64::INFINITY } else { f64::NEG_INFINITY },
                if name == "min" { f64::min } else { f64::max },
            ))
        }
        "clamp" => {
            if a.len() != 3 {
                return Err("clamp takes 3 arguments".into());
            }
            Ok(a[0].clamp(a[1], a[2]))
        }
        _ => Err(format!("unknown function {name:?}")),
    }
}
