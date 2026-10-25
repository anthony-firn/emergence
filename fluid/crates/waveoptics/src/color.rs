//! Colorimetry and image files.
//!
//! * Cone fundamentals: Stockman & Sharpe (Vision Res. 40, 1711 (2000)) 2°
//!   L, M, S in linear energy units, peak-normalised, 390–830 nm at 5 nm
//!   (CVRL file `linss2_10e_5.csv`).
//! * Colour matching: CIE (2006/2012) 2° XYZ functions defined from those
//!   cone fundamentals (CVRL `lin2012xyz2e_5_7sf.csv`). The LMS → XYZ
//!   matrix is recovered from the two tables by least squares (they are
//!   related by an exact linear transform, so the fit residual is the
//!   tables' rounding).
//! * Display: XYZ → linear sRGB (IEC 61966-2-1 matrix) and the sRGB
//!   transfer curve.

use std::io::BufWriter;
use std::path::Path;
use std::sync::OnceLock;

const CONES: &str = include_str!("../data/linss2_10e_5.csv");
const XYZ: &str = include_str!("../data/lin2012xyz2e_5_7sf.csv");

fn parse(text: &str) -> Vec<(f64, [f64; 3])> {
    text.lines()
        .filter_map(|l| {
            let f: Vec<&str> = l.split(',').map(|s| s.trim()).collect();
            let nm: f64 = f.first()?.parse().ok()?;
            let v = |i: usize| f.get(i).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
            Some((nm, [v(1), v(2), v(3)]))
        })
        .collect()
}

fn table(which: usize) -> &'static Vec<(f64, [f64; 3])> {
    static T: OnceLock<[Vec<(f64, [f64; 3])>; 2]> = OnceLock::new();
    &T.get_or_init(|| [parse(CONES), parse(XYZ)])[which]
}

fn interp(t: &[(f64, [f64; 3])], lambda: f64) -> [f64; 3] {
    let nm = lambda * 1e9;
    if nm < t[0].0 || nm > t[t.len() - 1].0 {
        return [0.0; 3];
    }
    let i = t.partition_point(|r| r.0 < nm).max(1);
    let (a, b) = (t[i - 1], t[i]);
    let f = (nm - a.0) / (b.0 - a.0);
    std::array::from_fn(|k| a.1[k] * (1.0 - f) + b.1[k] * f)
}

/// Cone fundamentals `(L, M, S)` at wavelength `λ` (m), energy units.
pub fn lms(lambda: f64) -> [f64; 3] {
    interp(table(0), lambda)
}

/// CIE 2006/2012 2° colour matching functions at `λ` (m).
pub fn xyz(lambda: f64) -> [f64; 3] {
    interp(table(1), lambda)
}

/// The LMS → XYZ matrix and the worst table residual of the fit.
pub fn lms_to_xyz() -> ([[f64; 3]; 3], f64) {
    static M: OnceLock<([[f64; 3]; 3], f64)> = OnceLock::new();
    *M.get_or_init(|| {
        let (c, x) = (table(0), table(1));
        let mut m = [[0.0; 3]; 3];
        // normal equations Aᵀ A, A rows = lms
        let mut ata = [[0.0; 3]; 3];
        let mut atb = [[0.0; 3]; 3];
        for (rc, rx) in c.iter().zip(x) {
            for i in 0..3 {
                for j in 0..3 {
                    ata[i][j] += rc.1[i] * rc.1[j];
                    atb[i][j] += rc.1[i] * rx.1[j];
                }
            }
        }
        let inv = invert3(ata);
        for out in 0..3 {
            for i in 0..3 {
                m[out][i] = (0..3).map(|k| inv[i][k] * atb[k][out]).sum();
            }
        }
        let mut worst: f64 = 0.0;
        for (rc, rx) in c.iter().zip(x) {
            for out in 0..3 {
                let p: f64 = (0..3).map(|i| m[out][i] * rc.1[i]).sum();
                worst = worst.max((p - rx.1[out]).abs());
            }
        }
        (m, worst)
    })
}

fn invert3(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    let mut inv = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (r1, r2) = ((j + 1) % 3, (j + 2) % 3);
            let (c1, c2) = ((i + 1) % 3, (i + 2) % 3);
            inv[i][j] = (a[r1][c1] * a[r2][c2] - a[r1][c2] * a[r2][c1]) / det;
        }
    }
    inv
}

/// XYZ → linear sRGB (IEC 61966-2-1, D65).
pub fn xyz_to_linear_srgb(v: [f64; 3]) -> [f64; 3] {
    [
        3.2406 * v[0] - 1.5372 * v[1] - 0.4986 * v[2],
        -0.9689 * v[0] + 1.8758 * v[1] + 0.0415 * v[2],
        0.0557 * v[0] - 0.2040 * v[1] + 1.0570 * v[2],
    ]
}

/// sRGB transfer function (linear 0–1 → encoded 0–1).
pub fn srgb_encode(c: f64) -> f64 {
    let c = c.clamp(0.0, 1.0);
    if c <= 0.003_130_8 {
        12.92 * c
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    }
}

/// Write an 8-bit sRGB PNG from linear RGB (scaled so `white` maps to 1),
/// out-of-gamut colours desaturated towards their luminance.
pub fn write_srgb_png(path: impl AsRef<Path>, width: usize, height: usize, rgb: &[[f64; 3]], white: f64) -> std::io::Result<()> {
    let mut data = Vec::with_capacity(width * height * 3);
    for p in rgb {
        let mut c = p.map(|v| v / white);
        let y = 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2];
        let lo = c.iter().cloned().fold(f64::INFINITY, f64::min);
        if lo < 0.0 && y > 0.0 {
            // desaturate until the most negative channel reaches zero
            let t = y / (y - lo);
            c = c.map(|v| y + t * (v - y));
        }
        for v in c {
            data.push((srgb_encode(v) * 255.0).round() as u8);
        }
    }
    let file = std::fs::File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), width as u32, height as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.set_source_srgb(png::SrgbRenderingIntent::Perceptual);
    let mut w = enc.write_header().map_err(std::io::Error::other)?;
    w.write_image_data(&data).map_err(std::io::Error::other)?;
    Ok(())
}

/// Write a 16-bit greyscale PNG of raw sensor values (no tone curve).
pub fn write_gray16_png(path: impl AsRef<Path>, width: usize, height: usize, v: &[u32], max: u32) -> std::io::Result<()> {
    let mut data = Vec::with_capacity(width * height * 2);
    for &x in v {
        let s = ((x as f64 / max.max(1) as f64) * 65535.0).round().clamp(0.0, 65535.0) as u16;
        data.extend_from_slice(&s.to_be_bytes());
    }
    let file = std::fs::File::create(path)?;
    let mut enc = png::Encoder::new(BufWriter::new(file), width as u32, height as u32);
    enc.set_color(png::ColorType::Grayscale);
    enc.set_depth(png::BitDepth::Sixteen);
    let mut w = enc.write_header().map_err(std::io::Error::other)?;
    w.write_image_data(&data).map_err(std::io::Error::other)?;
    Ok(())
}
