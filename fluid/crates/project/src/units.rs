//! Physical dimensions and units.
//!
//! Every number that crosses from one discipline to another carries its
//! dimension, so "0.5" can never be read as millimetres by one team and
//! inches by another. Values are stored in SI; units only matter when text
//! is parsed or printed. Parsing is strict: an unknown symbol, a mismatched
//! dimension or an affine unit (`degC`) inside a compound is an error, never
//! a guess.

use std::fmt;

/// SI base units, in the order of [`Dim`]'s exponents.
pub const BASE: [&str; 7] = ["m", "kg", "s", "A", "K", "mol", "cd"];

/// Exponents of the seven SI base units.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default, PartialOrd, Ord)]
pub struct Dim(pub [i8; 7]);

impl Dim {
    pub const NONE: Dim = Dim([0; 7]);
    pub const LENGTH: Dim = Dim([1, 0, 0, 0, 0, 0, 0]);
    pub const MASS: Dim = Dim([0, 1, 0, 0, 0, 0, 0]);
    pub const TIME: Dim = Dim([0, 0, 1, 0, 0, 0, 0]);
    pub const CURRENT: Dim = Dim([0, 0, 0, 1, 0, 0, 0]);
    pub const TEMPERATURE: Dim = Dim([0, 0, 0, 0, 1, 0, 0]);

    pub fn mul(self, o: Dim) -> Dim {
        let mut d = self.0;
        for (a, b) in d.iter_mut().zip(o.0) {
            *a += b;
        }
        Dim(d)
    }
    pub fn div(self, o: Dim) -> Dim {
        let mut d = self.0;
        for (a, b) in d.iter_mut().zip(o.0) {
            *a -= b;
        }
        Dim(d)
    }
    pub fn powi(self, n: i32) -> Dim {
        Dim(self.0.map(|a| (a as i32 * n) as i8))
    }
    /// Square root, if every exponent is even.
    pub fn sqrt(self) -> Option<Dim> {
        self.0.iter().all(|a| a % 2 == 0).then(|| Dim(self.0.map(|a| a / 2)))
    }
    pub fn is_none(self) -> bool {
        self == Dim::NONE
    }

    /// Plain SI spelling from base units, e.g. `kg·m^2·s^-3`.
    fn base_symbol(self) -> String {
        let mut parts = Vec::new();
        for (sym, e) in BASE.iter().zip(self.0) {
            match e {
                0 => {}
                1 => parts.push(sym.to_string()),
                e => parts.push(format!("{sym}^{e}")),
            }
        }
        if parts.is_empty() {
            "1".into()
        } else {
            parts.join("·")
        }
    }

    /// The most readable SI spelling: a named unit (`W`), a ratio or product
    /// of two or three named units (`K/W`, `W/(m·K)`), else base units.
    pub fn symbol(self) -> String {
        if self.is_none() {
            return String::new();
        }
        let named: Vec<(&str, Dim)> = NAMED.iter().map(|&(s, d)| (s, d)).collect();
        for &(s, d) in &named {
            if d == self {
                return s.into();
            }
        }
        for &(a, da) in &named {
            for &(b, db) in &named {
                if da.div(db) == self {
                    return format!("{a}/{b}");
                }
            }
        }
        for &(a, da) in &named {
            if da.powi(2) == self {
                return format!("{a}^2");
            }
            if da.powi(3) == self {
                return format!("{a}^3");
            }
            for &(b, db) in &named {
                if da.div(db.powi(2)) == self {
                    return format!("{a}/{b}^2");
                }
                for &(c, dc) in &named {
                    if da.div(db.mul(dc)) == self && b <= c {
                        return format!("{a}/({b}·{c})");
                    }
                }
            }
        }
        self.base_symbol()
    }
}

impl fmt::Display for Dim {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_none() {
            write!(f, "dimensionless")
        } else {
            write!(f, "{}", self.base_symbol())
        }
    }
}

const fn d(e: [i8; 7]) -> Dim {
    Dim(e)
}

/// Named SI units used for printing, most preferred first.
const NAMED: &[(&str, Dim)] = &[
    ("m", d([1, 0, 0, 0, 0, 0, 0])),
    ("kg", d([0, 1, 0, 0, 0, 0, 0])),
    ("s", d([0, 0, 1, 0, 0, 0, 0])),
    ("A", d([0, 0, 0, 1, 0, 0, 0])),
    ("K", d([0, 0, 0, 0, 1, 0, 0])),
    ("W", d([2, 1, -3, 0, 0, 0, 0])),
    ("J", d([2, 1, -2, 0, 0, 0, 0])),
    ("N", d([1, 1, -2, 0, 0, 0, 0])),
    ("Pa", d([-1, 1, -2, 0, 0, 0, 0])),
    ("V", d([2, 1, -3, -1, 0, 0, 0])),
    ("Ω", d([2, 1, -3, -2, 0, 0, 0])),
    ("Hz", d([0, 0, -1, 0, 0, 0, 0])),
    ("C", d([0, 0, 1, 1, 0, 0, 0])),
    ("F", d([-2, -1, 4, 2, 0, 0, 0])),
    ("S", d([-2, -1, 3, 2, 0, 0, 0])),
];

/// A unit: `SI value = value · scale + offset`.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Unit {
    pub scale: f64,
    pub offset: f64,
    pub dim: Dim,
}

impl Unit {
    pub const ONE: Unit = Unit { scale: 1.0, offset: 0.0, dim: Dim::NONE };

    pub fn to_si(&self, v: f64) -> f64 {
        v * self.scale + self.offset
    }
    pub fn from_si(&self, v: f64) -> f64 {
        (v - self.offset) / self.scale
    }
    /// True for units with an offset (`degC`), which only make sense for
    /// absolute values, never for differences or tolerances.
    pub fn is_affine(&self) -> bool {
        self.offset != 0.0
    }
}

/// Symbols recognised in unit strings: (symbol, scale to SI, dimension,
/// accepts an SI prefix). Exact matches win over prefixed ones, so `min`
/// is a minute and `mil` a thousandth of an inch, not milli-anything.
const SYMBOLS: &[(&str, f64, Dim, bool)] = &[
    ("m", 1.0, d([1, 0, 0, 0, 0, 0, 0]), true),
    ("g", 1e-3, d([0, 1, 0, 0, 0, 0, 0]), true),
    ("s", 1.0, d([0, 0, 1, 0, 0, 0, 0]), true),
    ("A", 1.0, d([0, 0, 0, 1, 0, 0, 0]), true),
    ("K", 1.0, d([0, 0, 0, 0, 1, 0, 0]), true),
    ("mol", 1.0, d([0, 0, 0, 0, 0, 1, 0]), true),
    ("Hz", 1.0, d([0, 0, -1, 0, 0, 0, 0]), true),
    ("N", 1.0, d([1, 1, -2, 0, 0, 0, 0]), true),
    ("Pa", 1.0, d([-1, 1, -2, 0, 0, 0, 0]), true),
    ("J", 1.0, d([2, 1, -2, 0, 0, 0, 0]), true),
    ("W", 1.0, d([2, 1, -3, 0, 0, 0, 0]), true),
    ("C", 1.0, d([0, 0, 1, 1, 0, 0, 0]), true),
    ("V", 1.0, d([2, 1, -3, -1, 0, 0, 0]), true),
    ("Ohm", 1.0, d([2, 1, -3, -2, 0, 0, 0]), true),
    ("Ω", 1.0, d([2, 1, -3, -2, 0, 0, 0]), true),
    ("S", 1.0, d([-2, -1, 3, 2, 0, 0, 0]), true),
    ("F", 1.0, d([-2, -1, 4, 2, 0, 0, 0]), true),
    ("L", 1e-3, d([3, 0, 0, 0, 0, 0, 0]), true),
    // Exact by definition (international inch, 1959).
    ("in", 0.0254, d([1, 0, 0, 0, 0, 0, 0]), false),
    ("mil", 2.54e-5, d([1, 0, 0, 0, 0, 0, 0]), false),
    ("ft", 0.3048, d([1, 0, 0, 0, 0, 0, 0]), false),
    ("lb", 0.45359237, d([0, 1, 0, 0, 0, 0, 0]), false),
    ("min", 60.0, d([0, 0, 1, 0, 0, 0, 0]), false),
    ("h", 3600.0, d([0, 0, 1, 0, 0, 0, 0]), false),
    ("%", 0.01, Dim::NONE, false),
    ("ppm", 1e-6, Dim::NONE, false),
];

const PREFIXES: &[(&str, f64)] = &[
    ("P", 1e15),
    ("T", 1e12),
    ("G", 1e9),
    ("M", 1e6),
    ("k", 1e3),
    ("c", 1e-2),
    ("m", 1e-3),
    ("u", 1e-6),
    ("µ", 1e-6),
    ("μ", 1e-6),
    ("n", 1e-9),
    ("p", 1e-12),
    ("f", 1e-15),
];

/// Absolute degrees Celsius: 0 degC = 273.15 K (exact, SI brochure).
const CELSIUS: Unit = Unit { scale: 1.0, offset: 273.15, dim: Dim::TEMPERATURE };

fn symbol(s: &str) -> Result<(f64, Dim), String> {
    if let Some(&(_, k, d, _)) = SYMBOLS.iter().find(|e| e.0 == s) {
        return Ok((k, d));
    }
    for &(p, f) in PREFIXES {
        if let Some(rest) = s.strip_prefix(p) {
            if let Some(&(_, k, d, _)) = SYMBOLS.iter().find(|e| e.0 == rest && e.3) {
                return Ok((f * k, d));
            }
        }
    }
    Err(format!("unknown unit `{s}`"))
}

/// Parse a unit expression: `mm`, `W/(m*K)`, `W/m/K`, `m^2`, `kg·m^-3`,
/// `degC` (alone only), `1` or the empty string (dimensionless).
pub fn parse_unit(text: &str) -> Result<Unit, String> {
    let t = text.trim();
    if t.is_empty() || t == "1" {
        return Ok(Unit::ONE);
    }
    if t == "degC" || t == "°C" {
        return Ok(CELSIUS);
    }
    if t.contains("degC") || t.contains("°C") {
        return Err(format!(
            "`{t}`: degC is an absolute temperature and is only valid alone; use K for differences"
        ));
    }
    let chars: Vec<char> = t.chars().collect();
    let mut p = UnitParser { c: &chars, i: 0 };
    let (k, d) = p.product()?;
    if p.i != chars.len() {
        return Err(format!("unit `{t}`: unexpected `{}`", chars[p.i..].iter().collect::<String>()));
    }
    Ok(Unit { scale: k, offset: 0.0, dim: d })
}

struct UnitParser<'a> {
    c: &'a [char],
    i: usize,
}

impl UnitParser<'_> {
    fn peek(&self) -> Option<char> {
        self.c.get(self.i).copied()
    }
    fn product(&mut self) -> Result<(f64, Dim), String> {
        let (mut k, mut d) = self.factor()?;
        loop {
            match self.peek() {
                Some('*') | Some('·') | Some(' ') => {
                    self.i += 1;
                    while self.peek() == Some(' ') {
                        self.i += 1;
                    }
                    if self.peek() == Some('/') {
                        continue;
                    }
                    let (k2, d2) = self.factor()?;
                    k *= k2;
                    d = d.mul(d2);
                }
                Some('/') => {
                    self.i += 1;
                    let (k2, d2) = self.factor()?;
                    k /= k2;
                    d = d.div(d2);
                }
                _ => return Ok((k, d)),
            }
        }
    }
    fn factor(&mut self) -> Result<(f64, Dim), String> {
        let (k, d) = if self.peek() == Some('(') {
            self.i += 1;
            let r = self.product()?;
            if self.peek() != Some(')') {
                return Err("unit: missing `)`".into());
            }
            self.i += 1;
            r
        } else {
            let start = self.i;
            while let Some(c) = self.peek() {
                if c.is_alphabetic() || c == '%' || c == 'µ' || c == 'Ω' {
                    self.i += 1;
                } else {
                    break;
                }
            }
            if start == self.i {
                if self.peek() == Some('1') {
                    self.i += 1;
                    (1.0, Dim::NONE)
                } else {
                    return Err(format!("unit: expected a symbol at `{}`", self.c[start..].iter().collect::<String>()));
                }
            } else {
                symbol(&self.c[start..self.i].iter().collect::<String>())?
            }
        };
        if self.peek() == Some('^') {
            self.i += 1;
            let start = self.i;
            if matches!(self.peek(), Some('-') | Some('+')) {
                self.i += 1;
            }
            while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                self.i += 1;
            }
            let e: i32 = self.c[start..self.i]
                .iter()
                .collect::<String>()
                .parse()
                .map_err(|_| "unit: bad exponent".to_string())?;
            return Ok((k.powi(e), d.powi(e)));
        }
        Ok((k, d))
    }
}

/// A number with a unit, already converted to SI.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Quantity {
    pub si: f64,
    pub dim: Dim,
}

/// Parse `"<number> <unit>"`, e.g. `"0.5 mm"`, `"40 degC"`, `"3 W/(m*K)"`.
/// Returns the quantity and the unit it was written in.
pub fn parse_quantity(text: &str) -> Result<(Quantity, Unit), String> {
    let (num, unit) = split_quantity(text)?;
    let unit = parse_unit(unit)?;
    Ok((Quantity { si: unit.to_si(num), dim: unit.dim }, unit))
}

/// Split `"<number> <unit>"` into the number and the unit text.
pub fn split_quantity(text: &str) -> Result<(f64, &str), String> {
    let t = text.trim();
    let end = t
        .char_indices()
        .find(|&(i, c)| {
            !(c.is_ascii_digit()
                || c == '.'
                || c == '+'
                || c == '-'
                || ((c == 'e' || c == 'E')
                    && t[i + 1..].starts_with(|n: char| n.is_ascii_digit() || n == '-' || n == '+')
                    && i > 0))
        })
        .map(|(i, _)| i)
        .unwrap_or(t.len());
    let num: f64 = t[..end].parse().map_err(|_| format!("`{t}`: expected a number first"))?;
    Ok((num, t[end..].trim()))
}

/// Format an SI value in `unit` (or plain SI) with `sig` significant digits.
pub fn format_si(si: f64, dim: Dim, unit: Option<&str>) -> String {
    match unit.map(parse_unit) {
        Some(Ok(u)) if u.dim == dim => {
            let s = unit.unwrap();
            format!("{} {s}", sig(u.from_si(si), 4))
        }
        _ => {
            let s = dim.symbol();
            if s.is_empty() {
                sig(si, 4)
            } else {
                format!("{} {s}", sig(si, 4))
            }
        }
    }
}

/// A number to `n` significant digits, without trailing zeros.
pub fn sig(v: f64, n: usize) -> String {
    if v == 0.0 || !v.is_finite() {
        return format!("{v}");
    }
    let mag = v.abs().log10().floor() as i32;
    if !(-4..6).contains(&mag) {
        let s = format!("{:.*e}", n.saturating_sub(1), v);
        // trim zeros in the mantissa
        if let Some((m, e)) = s.split_once('e') {
            let m = if m.contains('.') { m.trim_end_matches('0').trim_end_matches('.') } else { m };
            return format!("{m}e{e}");
        }
        return s;
    }
    let decimals = (n as i32 - 1 - mag).max(0) as usize;
    let s = format!("{v:.decimals$}");
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        s
    }
}
