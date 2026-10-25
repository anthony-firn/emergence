//! A synthesizable subset of Verilog-2001.
//!
//! Enough to write datapaths, controllers and small processors:
//!
//! ```verilog
//! module acc8(input clk, input rst, input [7:0] x, input [1:0] op,
//!             output [7:0] y, output reg [7:0] acc);
//!   wire [7:0] sum = acc + x;
//!   assign y = (op == 2'd0) ? sum : (op == 2'd1 ? acc & x : acc ^ x);
//!   always @(posedge clk) begin
//!     if (rst) acc <= 8'd0;
//!     else acc <= y;
//!   end
//! endmodule
//! ```
//!
//! Supported: one module per file with ANSI ports (`input`, `output`,
//! `output reg`, `[msb:lsb]`), `wire` and `reg` declarations (with an
//! optional initialiser on a wire), `assign`, `localparam`, and
//! `always @(posedge clk)` blocks of non-blocking assignments, `begin`/`end`,
//! `if`/`else` and `case` on constants. Expressions: sized and unsized
//! literals, bit and part selects, concatenation and replication, unary
//! `~ ! - & | ^`, binary `& | ^ + - * == != < <= > >= << >> && ||` (unsigned;
//! shifts and `*` by constants or small widths), and `?:`. Widths follow
//! Verilog's context rules in simplified form: an operator's result is as
//! wide as its widest operand (comparisons are one bit, `+` keeps the
//! widest operand's width unless the context is wider), and assignments
//! truncate or zero-extend.
//!
//! What is refused, with a message: latches (`always @*` with incomplete
//! assignment), multiple clocks, asynchronous resets, `initial`, tasks,
//! functions, generate blocks, signed arithmetic and anything else outside
//! the subset. A refusal is better than a silently different circuit.

use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    /// value, width (None: unsized)
    Const(u128, Option<usize>),
    Ident(String),
    /// signal, msb, lsb (bit select when msb == lsb)
    Slice(Box<Expr>, usize, usize),
    Concat(Vec<Expr>),
    Repeat(usize, Box<Expr>),
    Unary(&'static str, Box<Expr>),
    Binary(&'static str, Box<Expr>, Box<Expr>),
    Cond(Box<Expr>, Box<Expr>, Box<Expr>),
}

#[derive(Clone, Debug, PartialEq)]
pub enum Stmt {
    Assign(String, Option<(usize, usize)>, Expr),
    If(Expr, Vec<Stmt>, Vec<Stmt>),
    Case(Expr, Vec<(Vec<Expr>, Vec<Stmt>)>, Vec<Stmt>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Dir {
    Input,
    Output,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Port {
    pub name: String,
    pub dir: Dir,
    pub width: usize,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Module {
    pub name: String,
    pub ports: Vec<Port>,
    /// Every signal and its width (ports, wires, regs).
    pub widths: BTreeMap<String, usize>,
    /// Registers (assigned in the clocked block).
    pub regs: Vec<String>,
    /// Continuous assignments, in source order: (target, part, value).
    pub assigns: Vec<(String, Option<(usize, usize)>, Expr)>,
    /// The clock of the single `always @(posedge clk)` domain.
    pub clock: Option<String>,
    pub clocked: Vec<Stmt>,
}

#[derive(Clone, Debug, PartialEq)]
enum Tok {
    Id(String),
    Num(u128, Option<usize>),
    Sym(&'static str),
}

const SYMS: &[&str] = &[
    "<<", ">>", "<=", ">=", "==", "!=", "&&", "||", "(", ")", "[", "]", "{", "}", ";", ",", ":", "?", "=", "+", "-", "*", "&", "|", "^",
    "~", "!", "<", ">", "@", "#", ".",
];

fn lex(src: &str) -> Result<Vec<(Tok, usize)>, String> {
    let c: Vec<char> = src.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    let mut line = 1;
    while i < c.len() {
        let ch = c[i];
        if ch == '\n' {
            line += 1;
            i += 1;
        } else if ch.is_whitespace() {
            i += 1;
        } else if ch == '/' && c.get(i + 1) == Some(&'/') {
            while i < c.len() && c[i] != '\n' {
                i += 1;
            }
        } else if ch == '/' && c.get(i + 1) == Some(&'*') {
            i += 2;
            while i + 1 < c.len() && !(c[i] == '*' && c[i + 1] == '/') {
                if c[i] == '\n' {
                    line += 1;
                }
                i += 1;
            }
            i += 2;
        } else if ch.is_ascii_alphabetic() || ch == '_' || ch == '$' {
            let s = i;
            while i < c.len() && (c[i].is_ascii_alphanumeric() || c[i] == '_' || c[i] == '$') {
                i += 1;
            }
            out.push((Tok::Id(c[s..i].iter().collect()), line));
        } else if ch.is_ascii_digit() || ch == '\'' {
            // [size]'[base]digits, or plain decimal
            let s = i;
            while i < c.len() && (c[i].is_ascii_digit() || c[i] == '_') {
                i += 1;
            }
            let size_txt: String = c[s..i].iter().filter(|x| **x != '_').collect();
            if i < c.len() && c[i] == '\'' {
                i += 1;
                if i < c.len() && (c[i] == 's' || c[i] == 'S') {
                    return Err(format!("line {line}: signed literals are not supported"));
                }
                let base = c.get(i).copied().ok_or(format!("line {line}: literal without base"))?.to_ascii_lowercase();
                i += 1;
                let radix = match base {
                    'b' => 2,
                    'o' => 8,
                    'd' => 10,
                    'h' => 16,
                    _ => return Err(format!("line {line}: unknown base '{base}")),
                };
                let ds = i;
                while i < c.len() && (c[i].is_ascii_hexdigit() || c[i] == '_' || matches!(c[i], 'x' | 'X' | 'z' | 'Z')) {
                    i += 1;
                }
                let digits: String = c[ds..i].iter().filter(|x| **x != '_').collect();
                if digits.contains(['x', 'X', 'z', 'Z']) {
                    return Err(format!("line {line}: x/z literals are not synthesizable here"));
                }
                let v = u128::from_str_radix(&digits, radix).map_err(|_| format!("line {line}: bad literal {digits}"))?;
                let w = if size_txt.is_empty() { None } else { Some(size_txt.parse().map_err(|_| format!("line {line}: bad size"))?) };
                out.push((Tok::Num(v, w), line));
            } else {
                out.push((Tok::Num(size_txt.parse().map_err(|_| format!("line {line}: bad number"))?, None), line));
            }
        } else if let Some(sym) = SYMS.iter().find(|s| c[i..].iter().take(s.len()).copied().eq(s.chars())) {
            out.push((Tok::Sym(sym), line));
            i += sym.len();
        } else {
            return Err(format!("line {line}: unexpected character '{ch}'"));
        }
    }
    Ok(out)
}

struct P {
    t: Vec<(Tok, usize)>,
    i: usize,
    params: BTreeMap<String, u128>,
}

impl P {
    fn line(&self) -> usize {
        self.t.get(self.i).or(self.t.last()).map(|x| x.1).unwrap_or(0)
    }
    fn err<T>(&self, msg: impl std::fmt::Display) -> Result<T, String> {
        Err(format!("line {}: {msg}", self.line()))
    }
    fn peek(&self) -> Option<&Tok> {
        self.t.get(self.i).map(|x| &x.0)
    }
    fn is(&self, s: &str) -> bool {
        matches!(self.peek(), Some(Tok::Sym(x)) if *x == s) || matches!(self.peek(), Some(Tok::Id(x)) if x == s)
    }
    fn eat(&mut self, s: &str) -> bool {
        if self.is(s) {
            self.i += 1;
            true
        } else {
            false
        }
    }
    fn expect(&mut self, s: &str) -> Result<(), String> {
        if self.eat(s) {
            Ok(())
        } else {
            let found = format!("{:?}", self.peek());
            self.err(format!("expected `{s}`, found {found}"))
        }
    }
    fn ident(&mut self) -> Result<String, String> {
        match self.peek().cloned() {
            Some(Tok::Id(s)) => {
                self.i += 1;
                Ok(s)
            }
            other => self.err(format!("expected a name, found {other:?}")),
        }
    }
    fn const_expr(&mut self) -> Result<u128, String> {
        let e = self.expr()?;
        self.eval_const(&e)
    }
    fn eval_const(&self, e: &Expr) -> Result<u128, String> {
        Ok(match e {
            Expr::Const(v, _) => *v,
            Expr::Ident(n) => *self.params.get(n).ok_or(format!("line {}: `{n}` is not a constant", self.line()))?,
            Expr::Binary(op, a, b) => {
                let (a, b) = (self.eval_const(a)?, self.eval_const(b)?);
                match *op {
                    "+" => a + b,
                    "-" => a.wrapping_sub(b),
                    "*" => a * b,
                    "<<" => a << b,
                    ">>" => a >> b,
                    _ => return self.err(format!("operator `{op}` in a constant")),
                }
            }
            _ => return self.err("expected a constant"),
        })
    }
    /// `[msb:lsb]` → width, or 1.
    fn range(&mut self) -> Result<usize, String> {
        if self.eat("[") {
            let msb = self.const_expr()?;
            self.expect(":")?;
            let lsb = self.const_expr()?;
            self.expect("]")?;
            if lsb != 0 {
                return self.err("vectors must be declared [msb:0]");
            }
            Ok(msb as usize + 1)
        } else {
            Ok(1)
        }
    }

    // Expressions, by precedence (lowest first).
    fn expr(&mut self) -> Result<Expr, String> {
        let c = self.binary(0)?;
        if self.eat("?") {
            let a = self.expr()?;
            self.expect(":")?;
            let b = self.expr()?;
            return Ok(Expr::Cond(Box::new(c), Box::new(a), Box::new(b)));
        }
        Ok(c)
    }
    fn binary(&mut self, level: usize) -> Result<Expr, String> {
        const LEVELS: &[&[&str]] = &[
            &["||"],
            &["&&"],
            &["|"],
            &["^"],
            &["&"],
            &["==", "!="],
            &["<", "<=", ">", ">="],
            &["<<", ">>"],
            &["+", "-"],
            &["*"],
        ];
        if level == LEVELS.len() {
            return self.unary();
        }
        let mut e = self.binary(level + 1)?;
        loop {
            let op = LEVELS[level].iter().find(|o| self.is(o)).copied();
            match op {
                Some(op) => {
                    self.i += 1;
                    let r = self.binary(level + 1)?;
                    e = Expr::Binary(op, Box::new(e), Box::new(r));
                }
                None => return Ok(e),
            }
        }
    }
    fn unary(&mut self) -> Result<Expr, String> {
        for op in ["~", "!", "-", "&", "|", "^"] {
            if self.is(op) {
                self.i += 1;
                let e = self.unary()?;
                return Ok(Expr::Unary(op, Box::new(e)));
            }
        }
        self.postfix()
    }
    fn postfix(&mut self) -> Result<Expr, String> {
        let mut e = self.primary()?;
        while self.eat("[") {
            let msb = self.const_expr()? as usize;
            let lsb = if self.eat(":") { self.const_expr()? as usize } else { msb };
            self.expect("]")?;
            if lsb > msb {
                return self.err("part selects must be [msb:lsb]");
            }
            e = Expr::Slice(Box::new(e), msb, lsb);
        }
        Ok(e)
    }
    fn primary(&mut self) -> Result<Expr, String> {
        match self.peek().cloned() {
            Some(Tok::Num(v, w)) => {
                self.i += 1;
                Ok(Expr::Const(v, w))
            }
            Some(Tok::Id(n)) => {
                self.i += 1;
                if let Some(&v) = self.params.get(&n) {
                    return Ok(Expr::Const(v, None));
                }
                Ok(Expr::Ident(n))
            }
            Some(Tok::Sym("(")) => {
                self.i += 1;
                let e = self.expr()?;
                self.expect(")")?;
                Ok(e)
            }
            Some(Tok::Sym("{")) => {
                self.i += 1;
                let first = self.expr()?;
                if self.eat("{") {
                    let n = self.eval_const(&first)? as usize;
                    let inner = self.concat_body()?;
                    self.expect("}")?;
                    self.expect("}")?;
                    return Ok(Expr::Repeat(n, Box::new(inner)));
                }
                let mut parts = vec![first];
                while self.eat(",") {
                    parts.push(self.expr()?);
                }
                self.expect("}")?;
                Ok(Expr::Concat(parts))
            }
            other => self.err(format!("expected an expression, found {other:?}")),
        }
    }
    fn concat_body(&mut self) -> Result<Expr, String> {
        let mut parts = vec![self.expr()?];
        while self.eat(",") {
            parts.push(self.expr()?);
        }
        Ok(if parts.len() == 1 { parts.pop().unwrap() } else { Expr::Concat(parts) })
    }

    fn lvalue(&mut self) -> Result<(String, Option<(usize, usize)>), String> {
        let n = self.ident()?;
        if self.eat("[") {
            let msb = self.const_expr()? as usize;
            let lsb = if self.eat(":") { self.const_expr()? as usize } else { msb };
            self.expect("]")?;
            return Ok((n, Some((msb, lsb))));
        }
        Ok((n, None))
    }

    fn stmt(&mut self) -> Result<Vec<Stmt>, String> {
        if self.eat("begin") {
            let mut out = Vec::new();
            while !self.eat("end") {
                if self.peek().is_none() {
                    return self.err("missing `end`");
                }
                out.extend(self.stmt()?);
            }
            return Ok(out);
        }
        if self.eat("if") {
            self.expect("(")?;
            let c = self.expr()?;
            self.expect(")")?;
            let a = self.stmt()?;
            let b = if self.eat("else") { self.stmt()? } else { vec![] };
            return Ok(vec![Stmt::If(c, a, b)]);
        }
        if self.eat("case") {
            self.expect("(")?;
            let sel = self.expr()?;
            self.expect(")")?;
            let mut arms = Vec::new();
            let mut default = Vec::new();
            while !self.eat("endcase") {
                if self.eat("default") {
                    self.eat(":");
                    default = self.stmt()?;
                    continue;
                }
                let mut labels = vec![self.expr()?];
                while self.eat(",") {
                    labels.push(self.expr()?);
                }
                self.expect(":")?;
                arms.push((labels, self.stmt()?));
            }
            return Ok(vec![Stmt::Case(sel, arms, default)]);
        }
        let (n, part) = self.lvalue()?;
        if self.eat("=") {
            return self.err("use non-blocking `<=` in clocked blocks");
        }
        self.expect("<=")?;
        let e = self.expr()?;
        self.expect(";")?;
        Ok(vec![Stmt::Assign(n, part, e)])
    }
}

/// Parse one module.
pub fn parse(src: &str) -> Result<Module, String> {
    let mut p = P { t: lex(src)?, i: 0, params: BTreeMap::new() };
    let mut m = Module::default();
    p.expect("module")?;
    m.name = p.ident()?;
    if p.is("#") {
        return p.err("module parameters are not supported; use localparam");
    }
    p.expect("(")?;
    if !p.eat(")") {
        loop {
            let dir = if p.eat("input") {
                Dir::Input
            } else if p.eat("output") {
                Dir::Output
            } else {
                return p.err("ports must be declared ANSI-style: input/output");
            };
            let is_reg = p.eat("reg");
            p.eat("wire");
            let w = p.range()?;
            loop {
                let name = p.ident()?;
                m.ports.push(Port { name: name.clone(), dir, width: w });
                m.widths.insert(name.clone(), w);
                if is_reg {
                    m.regs.push(name);
                }
                if !(p.is(",") && !matches!(p.t.get(p.i + 1).map(|x| &x.0), Some(Tok::Id(s)) if s == "input" || s == "output")) {
                    break;
                }
                p.i += 1;
            }
            if !p.eat(",") {
                break;
            }
        }
        p.expect(")")?;
    }
    p.expect(";")?;
    while !p.eat("endmodule") {
        if p.peek().is_none() {
            return p.err("missing `endmodule`");
        }
        if p.eat("localparam") || p.eat("parameter") {
            p.range()?;
            loop {
                let n = p.ident()?;
                p.expect("=")?;
                let v = p.const_expr()?;
                p.params.insert(n, v);
                if !p.eat(",") {
                    break;
                }
            }
            p.expect(";")?;
        } else if p.is("wire") || p.is("reg") {
            let is_reg = p.eat("reg");
            if !is_reg {
                p.expect("wire")?;
            }
            let w = p.range()?;
            loop {
                let n = p.ident()?;
                if m.widths.contains_key(&n) {
                    return p.err(format!("`{n}` is declared twice"));
                }
                m.widths.insert(n.clone(), w);
                if is_reg {
                    m.regs.push(n.clone());
                }
                if p.eat("=") {
                    if is_reg {
                        return p.err("reg initialisers are not synthesizable; reset the register instead");
                    }
                    let e = p.expr()?;
                    m.assigns.push((n, None, e));
                }
                if !p.eat(",") {
                    break;
                }
            }
            p.expect(";")?;
        } else if p.eat("assign") {
            if p.eat("{") {
                // {a, b[3:0], c} = e: a hidden wire as wide as the targets,
                // sliced into each of them.
                let mut targets = vec![p.lvalue()?];
                while p.eat(",") {
                    targets.push(p.lvalue()?);
                }
                p.expect("}")?;
                p.expect("=")?;
                let e = p.expr()?;
                p.expect(";")?;
                let mut widths = Vec::new();
                for (n, part) in &targets {
                    let w = match part {
                        Some((a, b)) => a - b + 1,
                        None => *m.widths.get(n).ok_or(format!("line {}: assignment to undeclared `{n}`", p.line()))?,
                    };
                    widths.push(w);
                }
                let total: usize = widths.iter().sum();
                let hidden = format!("$concat{}", m.assigns.len());
                m.widths.insert(hidden.clone(), total);
                m.assigns.push((hidden.clone(), None, e));
                let mut lsb = total;
                for ((n, part), w) in targets.into_iter().zip(widths) {
                    lsb -= w;
                    m.assigns.push((n, part, Expr::Slice(Box::new(Expr::Ident(hidden.clone())), lsb + w - 1, lsb)));
                }
                continue;
            }
            let (n, part) = p.lvalue()?;
            p.expect("=")?;
            let e = p.expr()?;
            p.expect(";")?;
            m.assigns.push((n, part, e));
        } else if p.eat("always") {
            p.expect("@")?;
            if p.eat("*") || (p.eat("(") && p.is("*")) {
                return p.err("combinational always blocks are not supported; use assign");
            }
            if !p.eat("posedge") {
                return p.err("only always @(posedge clk) is supported");
            }
            let clk = p.ident()?;
            if p.is("or") || p.is(",") {
                return p.err("asynchronous resets are not supported; reset synchronously inside the block");
            }
            p.expect(")")?;
            if let Some(c) = &m.clock {
                if *c != clk {
                    return p.err(format!("a second clock `{clk}`: one clock domain is supported"));
                }
            }
            m.clock = Some(clk);
            let body = p.stmt()?;
            m.clocked.extend(body);
        } else if p.is("initial") {
            return p.err("initial blocks are not synthesizable");
        } else {
            let found = format!("{:?}", p.peek());
            return p.err(format!("unsupported construct at {found}"));
        }
    }
    // Every clocked target must be a declared reg, and every reg clocked.
    let mut targets = Vec::new();
    fn collect(s: &[Stmt], out: &mut Vec<String>) {
        for st in s {
            match st {
                Stmt::Assign(n, _, _) => {
                    if !out.contains(n) {
                        out.push(n.clone());
                    }
                }
                Stmt::If(_, a, b) => {
                    collect(a, out);
                    collect(b, out);
                }
                Stmt::Case(_, arms, d) => {
                    for (_, b) in arms {
                        collect(b, out);
                    }
                    collect(d, out);
                }
            }
        }
    }
    collect(&m.clocked, &mut targets);
    for t in &targets {
        if !m.regs.contains(t) {
            return Err(format!("`{t}` is assigned in the clocked block but not declared reg"));
        }
    }
    for r in &m.regs {
        if !targets.contains(r) {
            return Err(format!("reg `{r}` is never assigned in a clocked block (a latch or a constant?)"));
        }
    }
    for (n, _, _) in &m.assigns {
        if !m.widths.contains_key(n) {
            return Err(format!("assignment to undeclared `{n}`"));
        }
        if m.regs.contains(n) {
            return Err(format!("`{n}` is a reg and cannot be assigned continuously"));
        }
    }
    Ok(m)
}
