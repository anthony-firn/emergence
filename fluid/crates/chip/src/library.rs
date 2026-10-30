//! The cell library at an operating corner (supply voltage, junction
//! temperature), characterised from the transistor level on demand and
//! memoised: every combinational cell on the NLDM grid with held-out
//! certificates, and the flip-flop.

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use semiconductor::gates::{DeviceModel, Ekv};
use semiconductor::seq::{self, FlopTiming};
use semiconductor::timing::{build_library, Library};

/// Characterisation grid (input slew 10–90 %, load), and the held-out
/// points that certify interpolation (`analysis/gates.md`).
pub const SLEWS: [f64; 4] = [15e-12, 40e-12, 120e-12, 350e-12];
pub const LOADS: [f64; 3] = [1e-15, 4e-15, 12e-15];
pub const HELD: [(f64, f64); 4] = [(55e-12, 2e-15), (55e-12, 8e-15), (180e-12, 2e-15), (180e-12, 8e-15)];
pub const CLOCK_SLEWS: [f64; 3] = [20e-12, 80e-12, 200e-12];

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Corner {
    /// Supply (V).
    pub vdd: f64,
    /// Junction temperature (K).
    pub temp: f64,
}

impl Corner {
    pub const NOMINAL: Corner = Corner { vdd: 1.2, temp: 300.0 };
}

pub struct CellLib {
    pub corner: Corner,
    pub comb: Library,
    pub dff: FlopTiming,
}

type Memo = Mutex<HashMap<(u64, u64), Arc<CellLib>>>;

fn memo() -> &'static Memo {
    static M: std::sync::OnceLock<Memo> = std::sync::OnceLock::new();
    M.get_or_init(Default::default)
}

/// Bump when characterisation code changes in a way the constants below do
/// not capture, so stale cached libraries are never used.
const CACHE_VERSION: u32 = 1;

/// Everything a characterised library depends on, hashed (FNV-1a: stable
/// across Rust versions, unlike the standard hasher).
fn fingerprint() -> u64 {
    use semiconductor::gates::{N_TEMPERATURE, P_TEMPERATURE, TEMPERATURE_RANGE};
    let (n, p) = (Ekv::n(), Ekv::p());
    let m = DeviceModel::new(&n, &p);
    let text = format!(
        "{CACHE_VERSION} {n:?} {p:?} {N_TEMPERATURE:?} {P_TEMPERATURE:?} {TEMPERATURE_RANGE:?} {} {} {} {} {SLEWS:?} {LOADS:?} {HELD:?} {CLOCK_SLEWS:?}",
        m.p_scale, m.c_g, m.c_d, m.c_s
    );
    text.bytes().fold(0xcbf29ce484222325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100000001b3))
}

fn cache_path(corner: Corner) -> Option<std::path::PathBuf> {
    let dir = std::env::var_os("FLUID_CACHE")
        .map(std::path::PathBuf::from)
        .or_else(|| std::env::var_os("XDG_CACHE_HOME").map(|d| std::path::PathBuf::from(d).join("fluid")))
        .or_else(|| std::env::var_os("HOME").map(|d| std::path::PathBuf::from(d).join(".cache/fluid")))?;
    Some(dir.join(format!("chip-lib-{:016x}-{:016x}-{:016x}.json", fingerprint(), corner.vdd.to_bits(), corner.temp.to_bits())))
}

fn load_cached(corner: Corner) -> Option<CellLib> {
    let text = std::fs::read_to_string(cache_path(corner)?).ok()?;
    let (comb, dff): (Library, FlopTiming) = serde_json::from_str(&text).ok()?;
    Some(CellLib { corner, comb, dff })
}

fn store(lib: &CellLib) {
    let Some(path) = cache_path(lib.corner) else { return };
    let Some(dir) = path.parent() else { return };
    if std::fs::create_dir_all(dir).is_err() {
        return;
    }
    let Ok(text) = serde_json::to_string(&(&lib.comb, &lib.dff)) else { return };
    let tmp = path.with_extension(format!("tmp{}", std::process::id()));
    if std::fs::write(&tmp, text).is_ok() {
        let _ = std::fs::rename(&tmp, &path);
    }
}

/// Install a library characterised elsewhere, as the browser does: it has no
/// disk cache and characterising cells there would take minutes.
pub fn preload(lib: CellLib) {
    let key = (lib.corner.vdd.to_bits(), lib.corner.temp.to_bits());
    memo().lock().unwrap().insert(key, Arc::new(lib));
}

/// The library at `corner`: from memory, else the disk cache, else
/// characterised (and cached).
pub fn at(corner: Corner) -> Result<Arc<CellLib>, String> {
    if !(0.6..=1.5).contains(&corner.vdd) {
        return Err(format!("supply {:.3} V is outside 0.6–1.5 V for this 1.2 V technology", corner.vdd));
    }
    let key = (corner.vdd.to_bits(), corner.temp.to_bits());
    if let Some(l) = memo().lock().unwrap().get(&key) {
        return Ok(l.clone());
    }
    if let Some(l) = load_cached(corner) {
        let l = Arc::new(l);
        memo().lock().unwrap().insert(key, l.clone());
        return Ok(l);
    }
    let n = Ekv::n().at(corner.temp, corner.vdd)?;
    let p = Ekv::p().at(corner.temp, corner.vdd)?;
    let mut model = DeviceModel::new(&n, &p);
    model.vdd = corner.vdd;
    let comb = build_library(&model, &SLEWS, &LOADS, &HELD);
    let dff = seq::characterise(&model, &CLOCK_SLEWS, &LOADS);
    let lib = CellLib { corner, comb, dff };
    store(&lib);
    let lib = Arc::new(lib);
    memo().lock().unwrap().insert(key, lib.clone());
    Ok(lib)
}
