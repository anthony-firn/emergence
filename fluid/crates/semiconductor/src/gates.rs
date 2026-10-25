//! Standard cells at the gate rung of the electronics ladder
//! (`analysis/emergence_ladder.md`): INV, NAND2, NOR2, XOR2 and a full adder,
//! each a transistor netlist simulated with the compact MOSFET current law.
//!
//! The transistor rung here is the same all-region EKV compact law used to
//! certify the table model in `compact.rs` (there it is the held-out
//! reference for the PCHIP table; here it is the device model itself).
//! One device model (`Ekv`) is provided for each polarity, calibrated so
//! that the on-current matches the scaled 150 nm technology of
//! `analysis/semiconductor.md`: I_on,n = 330 µA/µm, I_on,p = 230 µA/µm,
//! V_dd = 1.2 V. The pMOS is width-scaled by `p_scale`.
//!
//! Everything transistor-level runs through `Net`: a flat list of MOSFETs
//! plus a capacitance per node. A cell is a small template; `flatten`
//! places templates into one flat net for the whole-circuit reference.
//! The transient is backward Euler with a scattered analytic Jacobian
//! (partials of each device current with respect to its three terminal
//! voltages) and damped Newton, so the series stacks are solved implicitly.
//!
//! Capacitances are the bias-averaged lumped values of
//! `analysis/semiconductor.md` §7: gate 1.031 fF/µm, drain junction
//! ≈ 1.4 fF/µm, source ≈ 0.7 fF/µm. They are the same lumped model the
//! ring oscillator uses, and are the largest approximation of this rung
//! (a charge-conserving terminal-charge model would be the next fold).

use crate::circuit::Fet;
use crate::timing::{GateNet, TimingTable};

/// Supply voltage of the scaled technology (analysis/semiconductor.md §6).
pub const VDD: f64 = 1.2;
/// Thermal voltage at 300 K, k_B T / q (CODATA / handbook).
pub const VT: f64 = 0.02585;
/// Temperatures (K) the device law's temperature dependence was fitted
/// over, and (dV_th/dT in V/K, I_s exponent) for each polarity
/// (`examples/ekv_temperature.rs`).
pub const TEMPERATURE_RANGE: (f64, f64) = (250.0, 425.0);
pub const N_TEMPERATURE: (f64, f64) = (-5.86e-4, 0.98);
pub const P_TEMPERATURE: (f64, f64) = (-6.20e-4, 0.94);
/// Minimum conductance to ground used only in the DC operating point, to pin
/// floating internal nodes of a stack (standard SPICE `gmin`).
pub const GMIN: f64 = 1e-12;

/// All-region EKV compact MOSFET current law.
///
/// I = I_s [ln²(1 + e^{v_p/2}) − ln²(1 + e^{(v_p − v_ds)/2})],
/// v_p = (V_gs − V_th)/(n V_T). This is the law `compact.rs` tests its
/// PCHIP table against. `pol` is +1 for nMOS, −1 for pMOS; `id` follows
/// the same sign convention as `compact::Table::id` (drain→source current,
/// negative for a pMOS pulling up).
#[derive(Clone, Copy, Debug)]
pub struct Ekv {
    pub isat: f64,
    pub vth: f64,
    pub slope: f64,
    pub pol: f64,
    /// Thermal voltage k_B T / q at the device temperature.
    pub vt: f64,
    /// Supply voltage: the largest |V_gs|, |V_ds| the law is used at.
    pub vmax: f64,
}

impl Ekv {
    #[inline]
    fn phi(x: f64) -> f64 {
        (1.0 + (0.5 * x).exp()).ln().powi(2)
    }

    /// Calibrate `isat` so that |I_d(V_dd, V_dd)| = `ion`.
    pub fn new(ion: f64, vth: f64, slope: f64, pol: f64) -> Self {
        let vp = (VDD - vth) / (slope * VT);
        let denom = Self::phi(vp) - Self::phi(vp - VDD / VT);
        Ekv { isat: ion / denom, vth, slope, pol, vt: VT, vmax: VDD }
    }

    /// Representative nMOS of the scaled 150 nm technology.
    pub fn n() -> Self {
        Self::new(330e-6, 0.32, 1.3, 1.0)
    }

    /// Representative pMOS (electrostatic mirror, lower mobility).
    pub fn p() -> Self {
        Self::new(230e-6, 0.32, 1.3, -1.0)
    }

    /// The same device at temperature `t` (K) and supply `vdd`.
    ///
    /// Folded from this crate's drift–diffusion transistors
    /// (`examples/ekv_temperature.rs`, `analysis/chip.md`): the threshold
    /// falls linearly and I_s scales as (T/300)^a (mobility and V_T²
    /// together), fitted to the TCAD ratios I_on(T)/I_on(300) and
    /// I_off(T)/I_off(300). `Err` outside the calibrated range.
    pub fn at(&self, t: f64, vdd: f64) -> Result<Ekv, String> {
        let (lo, hi) = TEMPERATURE_RANGE;
        if !(lo..=hi).contains(&t) {
            return Err(format!("{t:.1} K is outside the {lo:.0}–{hi:.0} K range the device law was calibrated over"));
        }
        let (dvth, a) = if self.pol > 0.0 { N_TEMPERATURE } else { P_TEMPERATURE };
        Ok(Ekv {
            isat: self.isat * (t / 300.0).powf(a),
            vth: self.vth + dvth * (t - 300.0),
            slope: self.slope,
            pol: self.pol,
            vt: VT * t / 300.0,
            vmax: vdd,
        })
    }

    /// Magnitude law on the n-sense positive (V_gs, V_ds) domain.
    #[inline]
    pub fn mag(&self, vgs: f64, vds: f64) -> f64 {
        let vp = (vgs - self.vth) / (self.slope * self.vt);
        let vd = vds.clamp(0.0, self.vmax);
        self.isat * (Self::phi(vp) - Self::phi(vp - vd / self.vt))
    }
}

impl Fet for Ekv {
    fn id(&self, vgs: f64, vds: f64) -> f64 {
        let (g, d) = (self.pol * vgs, self.pol * vds);
        let (g, d, sg) = if d >= 0.0 { (g, d, 1.0) } else { (g - d, -d, -1.0) };
        self.pol * sg * self.mag(g.clamp(-2.0, self.vmax), d)
    }
}

/// A transistor netlist device: `p` selects the pMOS law; `scale` multiplies
/// the current (width scaling).
#[derive(Clone, Copy, Debug)]
pub struct Mos {
    pub p: bool,
    pub g: usize,
    pub d: usize,
    pub s: usize,
    pub scale: f64,
}

/// A flat transistor network. Nodes 0 and 1 are ground and V_dd in every
/// template built here; `driven` marks nodes held by an ideal source.
#[derive(Clone, Debug, Default)]
pub struct Net {
    pub n_nodes: usize,
    pub driven: Vec<bool>,
    pub cap: Vec<f64>,
    pub devs: Vec<Mos>,
}

/// The device model: a unit-width nMOS and pMOS law, the pMOS width ratio,
/// and the lumped capacitances.
pub struct DeviceModel<'a> {
    pub n: &'a dyn Fet,
    pub p: &'a dyn Fet,
    pub p_scale: f64,
    pub c_g: f64,
    pub c_d: f64,
    pub c_s: f64,
    /// Supply voltage the cells run at.
    pub vdd: f64,
}

impl<'a> DeviceModel<'a> {
    pub fn new(n: &'a dyn Fet, p: &'a dyn Fet) -> Self {
        // W_p/W_n ≈ I_on,n/I_on,p for matched drive; analysis §7 uses 1.56
        // for a centred switching threshold. 1.5 is the nominal.
        DeviceModel { n, p, p_scale: 1.5, c_g: 1.031e-15, c_d: 1.4e-15, c_s: 0.7e-15, vdd: VDD }
    }

    /// Drain→source current for one device at absolute node voltages.
    #[inline]
    pub fn current(&self, dev: &Mos, vg: f64, vd: f64, vs: f64) -> f64 {
        if dev.p {
            self.p.id(vg - vs, vd - vs) * self.p_scale * dev.scale
        } else {
            self.n.id(vg - vs, vd - vs) * dev.scale
        }
    }

    /// Gate capacitance of one device (width-scaled).
    #[inline]
    pub fn gate_cap(&self, dev: &Mos) -> f64 {
        if dev.p {
            self.p_scale * dev.scale * self.c_g
        } else {
            dev.scale * self.c_g
        }
    }

    /// Drain/source diffusion capacitances (width-scaled).
    #[inline]
    pub fn drain_cap(&self, dev: &Mos) -> f64 {
        if dev.p { self.p_scale * dev.scale * self.c_d } else { dev.scale * self.c_d }
    }
    #[inline]
    pub fn source_cap(&self, dev: &Mos) -> f64 {
        if dev.p { self.p_scale * dev.scale * self.c_s } else { dev.scale * self.c_s }
    }
}

/// Cell types of the standard-cell library.
#[derive(Clone, Copy, PartialEq, Eq, Debug, serde::Serialize, serde::Deserialize)]
pub enum CellKind {
    Inv,
    Nand2,
    Nor2,
    Xor2,
    FullAdder,
}

impl CellKind {
    pub const ALL: [CellKind; 5] = [CellKind::Inv, CellKind::Nand2, CellKind::Nor2, CellKind::Xor2, CellKind::FullAdder];

    pub fn name(&self) -> &'static str {
        match self {
            CellKind::Inv => "INV",
            CellKind::Nand2 => "NAND2",
            CellKind::Nor2 => "NOR2",
            CellKind::Xor2 => "XOR2",
            CellKind::FullAdder => "FA",
        }
    }

    pub fn n_inputs(&self) -> usize {
        match self {
            CellKind::Inv => 1,
            CellKind::Nand2 | CellKind::Nor2 | CellKind::Xor2 => 2,
            CellKind::FullAdder => 3,
        }
    }

    pub fn n_outputs(&self) -> usize {
        if *self == CellKind::FullAdder { 2 } else { 1 }
    }

    /// Boolean function. Inputs and outputs are in `ins` order.
    pub fn eval(&self, ins: &[bool]) -> Vec<bool> {
        match self {
            CellKind::Inv => vec![!ins[0]],
            CellKind::Nand2 => vec![!(ins[0] && ins[1])],
            CellKind::Nor2 => vec![!(ins[0] || ins[1])],
            CellKind::Xor2 => vec![ins[0] ^ ins[1]],
            CellKind::FullAdder => vec![
                ins[0] ^ ins[1] ^ ins[2],
                (ins[0] && ins[1]) || (ins[0] && ins[2]) || (ins[1] && ins[2]),
            ],
        }
    }

    /// A local transistor template: nodes 0 = gnd, 1 = vdd.
    pub fn template(&self) -> Template {
        match self {
            CellKind::Inv => {
                let devs = vec![
                    Mos { p: false, g: 2, d: 3, s: 0, scale: 1.0 },
                    Mos { p: true, g: 2, d: 3, s: 1, scale: 1.0 },
                ];
                Template { n_nodes: 4, inputs: vec![2], outputs: vec![3], devs }
            }
            CellKind::Nand2 => {
                let devs = vec![
                    Mos { p: false, g: 2, d: 4, s: 5, scale: 1.0 },
                    Mos { p: false, g: 3, d: 5, s: 0, scale: 1.0 },
                    Mos { p: true, g: 2, d: 4, s: 1, scale: 1.0 },
                    Mos { p: true, g: 3, d: 4, s: 1, scale: 1.0 },
                ];
                Template { n_nodes: 6, inputs: vec![2, 3], outputs: vec![4], devs }
            }
            CellKind::Nor2 => {
                let devs = vec![
                    Mos { p: true, g: 2, d: 4, s: 5, scale: 1.0 },
                    Mos { p: true, g: 3, d: 5, s: 1, scale: 1.0 },
                    Mos { p: false, g: 2, d: 4, s: 0, scale: 1.0 },
                    Mos { p: false, g: 3, d: 4, s: 0, scale: 1.0 },
                ];
                Template { n_nodes: 6, inputs: vec![2, 3], outputs: vec![4], devs }
            }
            CellKind::Xor2 => {
                let mut t = Template { n_nodes: 2, inputs: vec![], outputs: vec![], devs: vec![] };
                let a = t.node();
                let b = t.node();
                let n1 = t.node();
                add_nand2(&mut t, a, b, n1);
                let n2 = t.node();
                add_nand2(&mut t, a, n1, n2);
                let n3 = t.node();
                add_nand2(&mut t, b, n1, n3);
                let y = t.node();
                add_nand2(&mut t, n2, n3, y);
                t.inputs = vec![a, b];
                t.outputs = vec![y];
                t
            }
            CellKind::FullAdder => {
                let mut t = Template { n_nodes: 2, inputs: vec![], outputs: vec![], devs: vec![] };
                let a = t.node();
                let b = t.node();
                let c = t.node();
                let t1 = t.node();
                add_xor2(&mut t, a, b, t1);
                let sum = t.node();
                add_xor2(&mut t, t1, c, sum);
                let n1 = t.node();
                add_nand2(&mut t, a, b, n1);
                let n2 = t.node();
                add_nand2(&mut t, t1, c, n2);
                let cout = t.node();
                add_nand2(&mut t, n1, n2, cout);
                t.inputs = vec![a, b, c];
                t.outputs = vec![sum, cout];
                t
            }
        }
    }
}

/// A local transistor template with local node indices.
#[derive(Clone, Debug)]
pub struct Template {
    pub n_nodes: usize,
    pub inputs: Vec<usize>,
    pub outputs: Vec<usize>,
    pub devs: Vec<Mos>,
}

impl Template {
    fn node(&mut self) -> usize {
        let k = self.n_nodes;
        self.n_nodes += 1;
        k
    }
}

/// NAND2 (transistors gated by `a`, `b`) driving local `out`; adds the
/// series-stack internal node.
fn add_nand2(t: &mut Template, a: usize, b: usize, out: usize) {
    let mid = t.node();
    t.devs.push(Mos { p: false, g: a, d: out, s: mid, scale: 1.0 });
    t.devs.push(Mos { p: false, g: b, d: mid, s: 0, scale: 1.0 });
    t.devs.push(Mos { p: true, g: a, d: out, s: 1, scale: 1.0 });
    t.devs.push(Mos { p: true, g: b, d: out, s: 1, scale: 1.0 });
}

/// XOR2 from four NAND2 cells (the classic 4-NAND XOR).
fn add_xor2(t: &mut Template, a: usize, b: usize, out: usize) {
    let n1 = t.node();
    add_nand2(t, a, b, n1);
    let n2 = t.node();
    add_nand2(t, a, n1, n2);
    let n3 = t.node();
    add_nand2(t, b, n1, n3);
    add_nand2(t, n2, n3, out);
}

/// Build the flat node capacitance list of a cell template: diffusion caps at
/// every drain/source, gate caps on every non-driven gate node. Gate nodes
/// that are cell inputs stay uncapacitated (they are driven).
pub fn net_from_template(model: &DeviceModel, t: &Template) -> Net {
    let mut net = Net { n_nodes: t.n_nodes, driven: vec![false; t.n_nodes], cap: vec![0.0; t.n_nodes], devs: t.devs.clone() };
    net.driven[0] = true;
    net.driven[1] = true;
    for &i in &t.inputs {
        net.driven[i] = true;
    }
    for dev in &t.devs {
        net.cap[dev.d] += model.drain_cap(dev);
        net.cap[dev.s] += model.source_cap(dev);
        if !net.driven[dev.g] {
            net.cap[dev.g] += model.gate_cap(dev);
        }
    }
    net
}

pub fn cell_net(model: &DeviceModel, kind: CellKind) -> Net {
    net_from_template(model, &kind.template())
}

/// Per-input-pin timing of one output: rise and fall NLDM tables, plus the
/// supply energy delivered on a rise/fall transition (which includes the
/// short-circuit term the pure ½CV² estimate misses).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct PinTiming {
    pub rise: TimingTable,
    pub fall: TimingTable,
    pub e_rise: TimingTable,
    pub e_fall: TimingTable,
}

/// Characterised timing and energy of one cell.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct CellTiming {
    pub kind: CellKind,
    /// input capacitance per cell input (F)
    pub c_in: Vec<f64>,
    /// self (drain-diffusion) capacitance seen at each output (F)
    pub c_out: Vec<f64>,
    /// [output][input pin] → rise/fall timing tables
    pub outs: Vec<Vec<PinTiming>>,
    /// worst-case leakage current drawn from V_dd (A)
    pub i_leak: f64,
}

/// A waveform: times and node voltages at each time.
#[derive(Clone, Debug)]
pub struct Wave {
    pub t: Vec<f64>,
    pub v: Vec<Vec<f64>>,
}

/// DC operating point (backward-Euler with the capacitive term dropped) of
/// the unknown nodes, given the driven node voltages. Newton on F=0.
pub fn dc_solve(net: &Net, model: &DeviceModel, drive: &dyn Fn(usize) -> f64, guess: &[f64]) -> Vec<f64> {
    let mut v = guess.to_vec();
    for i in 0..net.n_nodes {
        if net.driven[i] {
            v[i] = drive(i);
        }
    }
    let unk: Vec<usize> = (0..net.n_nodes).filter(|&i| !net.driven[i]).collect();
    if unk.is_empty() {
        return v;
    }
    let idx: Vec<isize> = {
        let mut m = vec![-1isize; net.n_nodes];
        for (k, &i) in unk.iter().enumerate() {
            m[i] = k as isize;
        }
        m
    };
    let nu = unk.len();
    for &i in &unk {
        v[i] = v[i].clamp(-0.5, model.vdd + 0.5);
    }
    for _ in 0..60 {
        let (f, mut j) = residual_and_jac(net, model, &v, &unk, &idx);
        // solve J d = -f
        let d = match solve(&mut j, &f.iter().map(|x| -x).collect::<Vec<_>>(), nu) {
            Some(d) => d,
            None => break,
        };
        let mut maxd = 0.0f64;
        for (k, &i) in unk.iter().enumerate() {
            let step = d[k].clamp(-0.5, 0.5);
            v[i] += step;
            maxd = maxd.max(step.abs());
        }
        if maxd < 1e-9 {
            break;
        }
    }
    v
}

/// Backward-Euler transient on the supplied (possibly non-uniform) time grid.
/// `drive(node, t)` gives the voltage of driven nodes. `extra_cap` adds a
/// capacitance at a node (the load). Returns the waveform.
pub fn transient(
    net: &Net,
    model: &DeviceModel,
    drive: &dyn Fn(usize, f64) -> f64,
    times: &[f64],
    v0: &[f64],
    extra_cap: &[(usize, f64)],
) -> Wave {
    let unk: Vec<usize> = (0..net.n_nodes).filter(|&i| !net.driven[i]).collect();
    let idx: Vec<isize> = {
        let mut m = vec![-1isize; net.n_nodes];
        for (k, &i) in unk.iter().enumerate() {
            m[i] = k as isize;
        }
        m
    };
    let nu = unk.len();
    let mut cap: Vec<f64> = net.cap.clone();
    for &(i, c) in extra_cap {
        cap[i] += c;
    }
    let mut v = v0.to_vec();
    for i in 0..net.n_nodes {
        if net.driven[i] {
            v[i] = drive(i, times[0]);
        }
    }
    let mut out = Wave { t: vec![times[0]], v: vec![v.clone()] };
    for w in times.windows(2) {
        let (t0, t1) = (w[0], w[1]);
        let dt = t1 - t0;
        let vold = v.clone();
        let mut ok = false;
        for _newton in 0..30 {
            let (f, mut j) = residual_and_jac_t(net, model, &v, &vold, drive, t1, dt, &unk, &idx, &cap);
            let d = match solve(&mut j, &f.iter().map(|x| -x).collect::<Vec<_>>(), nu) {
                Some(d) => d,
                None => break,
            };
            let mut maxd = 0.0f64;
            for (k, &i) in unk.iter().enumerate() {
                let step = d[k].clamp(-0.5, 0.5);
                v[i] += step;
                maxd = maxd.max(step.abs());
            }
            if maxd < 1e-7 {
                ok = true;
                break;
            }
        }
        // If Newton failed, fall back to explicit prediction (should be rare).
        if !ok {
            for &i in &unk {
                v[i] = vold[i];
            }
        }
        // store the correct driven-node voltages for energy accounting
        for i in 0..net.n_nodes {
            if net.driven[i] {
                v[i] = drive(i, t1);
            }
        }
        out.t.push(t1);
        out.v.push(v.clone());
    }
    out
}

/// Residual F and Jacobian dR/dV for the backward-Euler step. R_i =
/// cap_i (V_i − V_i,old)/dt − F_i(V). Only unknown rows/cols are filled.
#[allow(clippy::too_many_arguments)]
fn residual_and_jac_t(
    net: &Net,
    model: &DeviceModel,
    v: &[f64],
    vold: &[f64],
    drive: &dyn Fn(usize, f64) -> f64,
    t: f64,
    dt: f64,
    unk: &[usize],
    idx: &[isize],
    cap: &[f64],
) -> (Vec<f64>, Vec<f64>) {
    let nu = unk.len();
    let mut f = vec![0.0; nu];
    let mut j = vec![0.0; nu * nu];
    // driven node voltages
    let mut vv = v.to_vec();
    for i in 0..net.n_nodes {
        if net.driven[i] {
            vv[i] = drive(i, t);
        }
    }
    // F and dF/dV for the devices, then R = cap term − F, J = cap/dt − dF/dV.
    add_device_fj(net, model, &vv, idx, &mut f, &mut j, nu);
    for (k, &i) in unk.iter().enumerate() {
        f[k] = cap[i] * (vv[i] - vold[i]) / dt - f[k];
    }
    // J currently holds +dF/dV; the residual wants cap/dt − dF/dV.
    for x in j.iter_mut() {
        *x = -*x;
    }
    for (k, &i) in unk.iter().enumerate() {
        j[k * nu + k] += cap[i] / dt;
    }
    (f, j)
}

/// Residual F and Jacobian dF/dV for the DC case (cap term absent). A small
/// conductance to ground (`GMIN`) pins floating internal nodes.
fn residual_and_jac(net: &Net, model: &DeviceModel, v: &[f64], unk: &[usize], idx: &[isize]) -> (Vec<f64>, Vec<f64>) {
    let nu = unk.len();
    let mut f = vec![0.0; nu];
    let mut j = vec![0.0; nu * nu];
    add_device_fj(net, model, v, idx, &mut f, &mut j, nu);
    for (k, &i) in unk.iter().enumerate() {
        f[k] += GMIN * v[i];
        j[k * nu + k] += GMIN;
    }
    (f, j)
}

/// Scatter each device's current and terminal-voltage partials into F and J.
/// F_i = Σ (current into node i). R = cap term − F, so J = cap/dt·I − dF/dV.
fn add_device_fj(
    net: &Net,
    model: &DeviceModel,
    v: &[f64],
    idx: &[isize],
    f: &mut [f64],
    j: &mut [f64],
    nu: usize,
) {
    let h = 1e-4;
    for dev in &net.devs {
        let (vg, vd, vs) = (v[dev.g], v[dev.d], v[dev.s]);
        let i0 = model.current(dev, vg, vd, vs);
        // Partial of the device current with respect to V_gs, V_ds.
        let (vgs, vds) = (vg - vs, vd - vs);
        let dvg = (model.id_dev(dev, vgs + h, vds) - model.id_dev(dev, vgs - h, vds)) / (2.0 * h);
        let dvd = (model.id_dev(dev, vgs, vds + h) - model.id_dev(dev, vgs, vds - h)) / (2.0 * h);
        // dI/dVs = -(dI/dVgs + dI/dVds)
        let dvs = -(dvg + dvd);
        let terms = [
            (dev.d, dev.g, -dvg),
            (dev.d, dev.d, -dvd),
            (dev.d, dev.s, -dvs),
            (dev.s, dev.g, dvg),
            (dev.s, dev.d, dvd),
            (dev.s, dev.s, dvs),
        ];
        let fd = idx[dev.d];
        let fs = idx[dev.s];
        if fd >= 0 {
            f[fd as usize] -= i0;
        }
        if fs >= 0 {
            f[fs as usize] += i0;
        }
        for &(r, c, val) in &terms {
            let (ri, ci) = (idx[r], idx[c]);
            if ri >= 0 && ci >= 0 {
                j[ri as usize * nu + ci as usize] += val;
            }
        }
    }
}

impl DeviceModel<'_> {
    #[inline]
    fn id_dev(&self, dev: &Mos, vgs: f64, vds: f64) -> f64 {
        if dev.p {
            self.p.id(vgs, vds) * self.p_scale * dev.scale
        } else {
            self.n.id(vgs, vds) * dev.scale
        }
    }
}

/// Solve a dense linear system A x = b by LU with partial pivoting.
/// `a` is overwritten by its LU factors.
pub fn solve(a: &mut [f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    // factor
    for k in 0..n {
        let mut piv = k;
        let mut mx = a[k * n + k].abs();
        for i in k + 1..n {
            let x = a[i * n + k].abs();
            if x > mx {
                mx = x;
                piv = i;
            }
        }
        if mx < 1e-300 {
            return None;
        }
        if piv != k {
            for col in 0..n {
                a.swap(k * n + col, piv * n + col);
            }
        }
        let akk = a[k * n + k];
        for i in k + 1..n {
            let f = a[i * n + k] / akk;
            a[i * n + k] = f;
            for col in k + 1..n {
                a[i * n + col] -= f * a[k * n + col];
            }
        }
    }
    let mut x = b.to_vec();
    for i in 0..n {
        let mut s = x[i];
        for k in 0..i {
            s -= a[i * n + k] * x[k];
        }
        x[i] = s;
    }
    for i in (0..n).rev() {
        let mut s = x[i];
        for k in i + 1..n {
            s -= a[i * n + k] * x[k];
        }
        x[i] = s / a[i * n + i];
    }
    Some(x)
}

// ---------------------------------------------------------------------------
// Transistor-level characterisation of the standard cells
// ---------------------------------------------------------------------------

/// One transistor-level characterisation result.
#[derive(Clone, Copy, Debug)]
pub struct CharPoint {
    /// propagation delay, input 50 % to output 50 % (s)
    pub delay: f64,
    /// output slew, 10–90 % rise or 90–10 % fall (s)
    pub out_slew: f64,
    /// energy delivered by V_dd over the transition (J)
    pub e_supply: f64,
}

/// Find a static setting of the other inputs and the initial value `p0` of
/// the switching pin `pin` such that flipping `pin` flips output `out` in
/// the requested direction.
fn static_for(kind: CellKind, pin: usize, out: usize, want_rise: bool) -> (Vec<bool>, bool) {
    let k = kind.n_inputs();
    for mask in 0..(1u32 << (k - 1)) {
        let mut others = vec![false; k];
        let mut m = 0;
        for (i, o) in others.iter_mut().enumerate() {
            if i != pin {
                *o = (mask >> m) & 1 == 1;
                m += 1;
            }
        }
        for p0 in [false, true] {
            let mut i0 = others.clone();
            i0[pin] = p0;
            let mut i1 = others.clone();
            i1[pin] = !p0;
            let (o0, o1) = (kind.eval(&i0)[out], kind.eval(&i1)[out]);
            let want = if want_rise { (false, true) } else { (true, false) };
            if (o0, o1) == want {
                return (others, p0);
            }
        }
    }
    panic!("no static setting flips output {out} from pin {pin}");
}

fn time_grid(t_start: f64, tr: f64, t_end: f64) -> Vec<f64> {
    let mut t = vec![t_start];
    let push = |t: &mut Vec<f64>, target: f64, dt: f64| {
        while *t.last().unwrap() < target - 1e-18 {
            let nt = (*t.last().unwrap() + dt).min(target);
            t.push(nt);
        }
    };
    let dt_e = (tr / 10.0).clamp(0.05e-12, 1e-12);
    push(&mut t, t_start + tr, dt_e);
    push(&mut t, t_end, 2e-12);
    t
}

fn crossing(w: &Wave, node: usize, level: f64, rising: bool, after: f64) -> f64 {
    for k in 1..w.t.len() {
        if w.t[k] < after {
            continue;
        }
        let (a, b) = (w.v[k - 1][node], w.v[k][node]);
        let hit = if rising { a < level && b >= level } else { a > level && b <= level };
        if hit {
            return w.t[k - 1] + (level - a) / (b - a) * (w.t[k] - w.t[k - 1]);
        }
    }
    f64::NAN
}

/// Current delivered by the V_dd source into the network at a node-voltage
/// vector (A).
pub fn supply_current(net: &Net, model: &DeviceModel, v: &[f64]) -> f64 {
    let mut i = 0.0;
    for dev in &net.devs {
        let cur = model.current(dev, v[dev.g], v[dev.d], v[dev.s]);
        if dev.s == 1 {
            i -= cur;
        }
        if dev.d == 1 {
            i += cur;
        }
    }
    i
}

/// Characterise one (input pin, output, edge) at one (input slew, load).
pub fn char_point(model: &DeviceModel, kind: CellKind, pin: usize, out: usize, rise: bool, s_in: f64, cl: f64) -> CharPoint {
    let tmpl = kind.template();
    let net = net_from_template(model, &tmpl);
    let (others, p0) = static_for(kind, pin, out, rise);
    let tr = s_in / 0.8;
    let t_start = 2.0 * tr + 20e-12;
    let t_in50 = t_start + 0.5 * tr;
    // tail long enough to charge the load through a few RC times
    let t_tail = 4.0 * cl * model.vdd / 2.0e-4 + 100e-12;
    let t_end = t_start + tr + (4.0 * tr).max(120e-12) + t_tail;
    let times = time_grid(t_start, tr, t_end);
    let pin_node = tmpl.inputs[pin];
    let in_nodes: Vec<usize> = tmpl.inputs.clone();
    let others_c = others.clone();
    let a0 = if p0 { model.vdd } else { 0.0 };
    let a1 = if p0 { 0.0 } else { model.vdd };
    let ramp = move |t: f64| -> f64 {
        if t <= t_start {
            a0
        } else if t >= t_start + tr {
            a1
        } else {
            a0 + (a1 - a0) * (t - t_start) / tr
        }
    };
    let drive = move |node: usize, t: f64| -> f64 {
        if node == 0 {
            0.0
        } else if node == 1 {
            model.vdd
        } else if node == pin_node {
            ramp(t)
        } else if let Some(i) = in_nodes.iter().position(|&x| x == node) {
            if others_c[i] {
                model.vdd
            } else {
                0.0
            }
        } else {
            0.5 * model.vdd
        }
    };
    // Initial guess: output nodes at their logic value, others at V_dd/2.
    // The circuit is settled over [0, t_start] by the same backward-Euler
    // transient, so no separate (and fragile) DC solve is needed: the node
    // capacitances define the floating stack nodes.
    let mut iv0 = others.clone();
    iv0[pin] = p0;
    let ov0 = kind.eval(&iv0);
    let mut guess = vec![0.5 * model.vdd; net.n_nodes];
    for (o, &on) in tmpl.outputs.iter().enumerate() {
        guess[on] = if ov0[o] { model.vdd } else { 0.0 };
    }
    // Coarse pre-settle (backward Euler, large steps) so that all internal
    // nodes start near the static state before the edge; the node
    // capacitances define the floating stack nodes.
    let v0 = settle(&net, model, &|n| drive(n, 0.0), &guess);
    let out_node = tmpl.outputs[out];
    let wave = transient(&net, model, &drive, &times, &v0, &[(out_node, cl)]);
    let t50 = crossing(&wave, out_node, 0.5 * model.vdd, rise, t_start);
    let delay = t50 - t_in50;
    let out_slew = if rise {
        crossing(&wave, out_node, 0.9 * model.vdd, true, t_start) - crossing(&wave, out_node, 0.1 * model.vdd, true, t_start)
    } else {
        crossing(&wave, out_node, 0.1 * model.vdd, false, t_start) - crossing(&wave, out_node, 0.9 * model.vdd, false, t_start)
    };
    // energy delivered by the supply over the transition
    let mut q = 0.0;
    for k in 1..wave.t.len() {
        if wave.t[k - 1] < t_start || wave.t[k] > t_end {
            continue;
        }
        let i0 = supply_current(&net, model, &wave.v[k - 1]);
        let i1 = supply_current(&net, model, &wave.v[k]);
        q += 0.5 * (i0 + i1) * (wave.t[k] - wave.t[k - 1]);
    }
    CharPoint { delay, out_slew, e_supply: model.vdd * q }
}

/// Settle to the DC operating point through the capacitive network (a
/// backward-Euler march with moderate steps) and return the final node
/// voltages. Used before the characterisation edge, and for leakage, where a
/// direct DC Newton can be ill-conditioned on floating stack nodes.
pub fn settle(net: &Net, model: &DeviceModel, drive: &dyn Fn(usize) -> f64, guess: &[f64]) -> Vec<f64> {
    let times: Vec<f64> = (0..=100).map(|k| k as f64 * 5e-12).collect();
    let drv = |n: usize, _t: f64| drive(n);
    let w = transient(net, model, &drv, &times, guess, &[]);
    w.v.last().unwrap().clone()
}

/// Worst-case leakage current (A) over all input states.
///
/// A direct DC Newton is ill-conditioned on the floating mid-nodes of series
/// stacks: the EKV off-current is ~1e-10 A, below the round-off of the on
/// current, and Newton runs to unphysical voltages (tests showed ±30 V).
/// Instead each transistor is first treated as a switch (on: 10 mS, off:
/// 1e-14 S), the node voltages of that linear network are solved, switches
/// are re-classified from the solved gate–source voltages and the loop is
/// repeated to a fixed point. The supply current is then evaluated with the
/// full non-linear law at those node voltages. This is the standard
/// switch-level leakage estimate and keeps every node in [0, V_dd].
fn leakage(model: &DeviceModel, kind: CellKind) -> f64 {
    let tmpl = kind.template();
    let net = net_from_template(model, &tmpl);
    let ni = kind.n_inputs();
    let mut worst = 0.0f64;
    for mask in 0..(1u32 << ni) {
        let iv: Vec<bool> = (0..ni).map(|i| (mask >> i) & 1 == 1).collect();
        // Fixed (driven) node voltages: ground, supply and the primary inputs.
        let fixed = |node: usize| -> Option<f64> {
            if node == 0 {
                Some(0.0)
            } else if node == 1 {
                Some(model.vdd)
            } else {
                tmpl.inputs.iter().position(|&x| x == node).map(|p| if iv[p] { model.vdd } else { 0.0 })
            }
        };
        let unk: Vec<usize> = (0..net.n_nodes).filter(|&i| fixed(i).is_none()).collect();
        let nu = unk.len();
        let idx: Vec<isize> = {
            let mut m = vec![-1isize; net.n_nodes];
            for (q, &i) in unk.iter().enumerate() {
                m[i] = q as isize;
            }
            m
        };
        let mut v = vec![0.0f64; net.n_nodes];
        for i in 0..net.n_nodes {
            if let Some(x) = fixed(i) {
                v[i] = x;
            }
        }
        for _pass in 0..40 {
            let mut g = vec![0.0; nu * nu];
            let mut b = vec![0.0; nu];
            for d in &net.devs {
                let on = if d.p { v[d.s] - v[d.g] > 0.2 } else { v[d.g] - v[d.s] > 0.2 };
                let gg = if on { 1e-2 } else { 1e-14 };
                let (rd, rs) = (idx[d.d], idx[d.s]);
                if rd >= 0 {
                    g[rd as usize * nu + rd as usize] += gg;
                    if let Some(s) = fixed(d.s) {
                        b[rd as usize] += gg * s;
                    }
                }
                if rs >= 0 {
                    g[rs as usize * nu + rs as usize] += gg;
                    if let Some(x) = fixed(d.d) {
                        b[rs as usize] += gg * x;
                    }
                }
                if rd >= 0 && rs >= 0 {
                    g[rd as usize * nu + rs as usize] -= gg;
                    g[rs as usize * nu + rd as usize] -= gg;
                }
            }
            let vv = match solve(&mut g, &b, nu) {
                Some(x) => x,
                None => break,
            };
            let mut maxd = 0.0f64;
            for (q, &i) in unk.iter().enumerate() {
                maxd = maxd.max((vv[q] - v[i]).abs());
                v[i] = vv[q];
            }
            if maxd < 1e-13 {
                break;
            }
        }
        worst = worst.max(supply_current(&net, model, &v).abs());
    }
    worst
}

/// Characterise a cell into NLDM tables on the given (input slew, load) grid,
/// with input capacitances and worst-case leakage.
pub fn characterise(model: &DeviceModel, kind: CellKind, slews: &[f64], loads: &[f64]) -> CellTiming {
    let tmpl = kind.template();
    let ni = kind.n_inputs();
    let no = kind.n_outputs();
    let mut c_in = vec![0.0; ni];
    for dev in &tmpl.devs {
        if let Some(p) = tmpl.inputs.iter().position(|&x| x == dev.g) {
            c_in[p] += model.gate_cap(dev);
        }
    }
    let mut c_out = vec![0.0; no];
    for (o, &on) in tmpl.outputs.iter().enumerate() {
        for dev in &tmpl.devs {
            if dev.d == on {
                c_out[o] += model.drain_cap(dev);
            }
        }
    }
    let mut outs = Vec::new();
    for o in 0..no {
        let mut pins = Vec::new();
        for p in 0..ni {
            let make = |rise: bool| -> (TimingTable, TimingTable) {
                let mut delay = vec![vec![0.0; loads.len()]; slews.len()];
                let mut sls = vec![vec![0.0; loads.len()]; slews.len()];
                let mut en = vec![vec![0.0; loads.len()]; slews.len()];
                // Grid points are independent transient runs.
                use rayon::prelude::*;
                let pts: Vec<(usize, usize)> = (0..slews.len()).flat_map(|i| (0..loads.len()).map(move |j| (i, j))).collect();
                let cps: Vec<CharPoint> = pts.par_iter().map(|&(i, j)| char_point(model, kind, p, o, rise, slews[i], loads[j])).collect();
                for (&(i, j), cp) in pts.iter().zip(cps) {
                    delay[i][j] = cp.delay;
                    sls[i][j] = cp.out_slew;
                    en[i][j] = cp.e_supply;
                }
                (
                    TimingTable::new(slews.to_vec(), loads.to_vec(), delay, sls),
                    TimingTable::new(slews.to_vec(), loads.to_vec(), en, vec![vec![0.0; loads.len()]; slews.len()]),
                )
            };
            let (rise, e_rise) = make(true);
            let (fall, e_fall) = make(false);
            pins.push(PinTiming { rise, fall, e_rise, e_fall });
        }
        outs.push(pins);
    }
    CellTiming { kind, c_in, c_out, outs, i_leak: leakage(model, kind) }
}

/// Attach the interpolation certificate to every table of a characterised
/// cell, from transistor-level held-out points (`(slew, load)` pairs that
/// were not on the table grid).
pub fn set_certificates(model: &DeviceModel, ct: &mut CellTiming, held: &[(f64, f64)]) {
    let kind = ct.kind;
    let no = kind.n_outputs();
    let ni = kind.n_inputs();
    let mut worst_e = 0.0f64;
    let mut worst_s = 0.0f64;
    for o in 0..no {
        for p in 0..ni {
            for rise in [true, false] {
                use rayon::prelude::*;
                let pts: Vec<(f64, f64, f64, f64)> = held
                    .par_iter()
                    .map(|&(s, l)| {
                        let cp = char_point(model, kind, p, o, rise, s, l);
                        (s, l, cp.delay, cp.out_slew)
                    })
                    .collect();
                let tab = if rise { &mut ct.outs[o][p].rise } else { &mut ct.outs[o][p].fall };
                let (e, se) = tab.certify(&pts);
                tab.cert_error = e;
                tab.cert_slew = se;
                worst_e = worst_e.max(e);
                worst_s = worst_s.max(se);
            }
        }
    }
    let _ = (worst_e, worst_s);
}

/// Flatten a gate-level netlist into one transistor net. Returns the net and
/// the global node index of every net.
pub fn flatten(model: &DeviceModel, gnet: &GateNet) -> (Net, Vec<usize>) {
    let mut net = Net { n_nodes: 2, driven: vec![true, true], cap: vec![0.0, 0.0], devs: vec![] };
    let mut map = vec![0usize; gnet.n_nets];
    for n in 0..gnet.n_nets {
        map[n] = net.n_nodes;
        net.n_nodes += 1;
        net.driven.push(false);
        net.cap.push(0.0);
    }
    for &n in &gnet.inputs {
        net.driven[map[n]] = true;
    }
    for inst in &gnet.insts {
        let t = inst.kind.template();
        let mut lm = vec![usize::MAX; t.n_nodes];
        lm[0] = 0;
        lm[1] = 1;
        for (i, &n) in inst.ins.iter().enumerate() {
            lm[t.inputs[i]] = map[n];
        }
        for (o, &n) in inst.outs.iter().enumerate() {
            lm[t.outputs[o]] = map[n];
        }
        for k in 2..t.n_nodes {
            if lm[k] == usize::MAX {
                lm[k] = net.n_nodes;
                net.n_nodes += 1;
                net.driven.push(false);
                net.cap.push(0.0);
            }
        }
        for d in &t.devs {
            let dev = Mos { p: d.p, g: lm[d.g], d: lm[d.d], s: lm[d.s], scale: d.scale };
            net.cap[dev.d] += model.drain_cap(&dev);
            net.cap[dev.s] += model.source_cap(&dev);
            if !net.driven[dev.g] {
                net.cap[dev.g] += model.gate_cap(&dev);
            }
            net.devs.push(dev);
        }
    }
    (net, map)
}

/// A primary-input stimulus: a list of (50 % crossing time, input vector).
/// Between steps the inputs are held; at a step each changed input ramps
/// linearly over the 10–90 % slew `slew`.
#[derive(Clone, Debug)]
pub struct Stimulus {
    pub steps: Vec<(f64, Vec<bool>)>,
    pub slew: f64,
    pub vdd: f64,
}

impl Stimulus {
    pub fn new(steps: Vec<(f64, Vec<bool>)>, slew: f64) -> Self {
        Stimulus { steps, slew, vdd: VDD }
    }

    pub fn with_vdd(mut self, vdd: f64) -> Self {
        self.vdd = vdd;
        self
    }

    /// Full 0–100 % ramp time of a linear edge with this 10–90 % slew.
    pub fn full_ramp(&self) -> f64 {
        self.slew / 0.8
    }

    /// Voltage of primary input `pin` at time `t`.
    pub fn value(&self, pin: usize, t: f64) -> f64 {
        let full = self.full_ramp();
        if t <= self.steps[0].0 {
            return if self.steps[0].1[pin] { self.vdd } else { 0.0 };
        }
        for k in 1..self.steps.len() {
            let (tk, to) = (&self.steps[k].0, &self.steps[k].1);
            let from = &self.steps[k - 1].1;
            if t < tk - 0.5 * full {
                return if from[pin] { self.vdd } else { 0.0 };
            }
            if t <= tk + 0.5 * full {
                let frac = (t - (tk - 0.5 * full)) / full;
                let a = if from[pin] { self.vdd } else { 0.0 };
                let b = if to[pin] { self.vdd } else { 0.0 };
                return a + (b - a) * frac;
            }
        }
        let last = &self.steps.last().unwrap().1;
        if last[pin] { self.vdd } else { 0.0 }
    }
}

/// Public 50 % (or arbitrary level) crossing of a node waveform.
pub fn wave_crossing(w: &Wave, node: usize, level: f64, rising: bool, after: f64) -> f64 {
    crossing(w, node, level, rising, after)
}

/// Simulate a gate-level netlist at transistor level (flattened), driving the
/// primary inputs with `stim`. Returns the flat waveform, the global node of
/// every net, and the flat net (for energy accounting).
pub struct FlatSim {
    pub wave: Wave,
    pub map: Vec<usize>,
    pub net: Net,
}

pub fn simulate_flat(model: &DeviceModel, gnet: &GateNet, stim: &Stimulus, times: &[f64]) -> FlatSim {
    let (fnet, map) = flatten(model, gnet);
    let all = gnet.eval(&stim.steps[0].1);
    let mut guess = vec![0.5 * model.vdd; fnet.n_nodes];
    for nn in 0..gnet.n_nets {
        guess[map[nn]] = if all[nn] { model.vdd } else { 0.0 };
    }
    let drive0 = |node: usize| -> f64 {
        if node == 0 {
            0.0
        } else if node == 1 {
            model.vdd
        } else {
            for (i, &n) in gnet.inputs.iter().enumerate() {
                if map[n] == node {
                    return if stim.steps[0].1[i] { model.vdd } else { 0.0 };
                }
            }
            0.5 * model.vdd
        }
    };
    let v0 = settle(&fnet, model, &drive0, &guess);
    let drive = |node: usize, t: f64| -> f64 {
        if node == 0 {
            0.0
        } else if node == 1 {
            model.vdd
        } else {
            for (i, &n) in gnet.inputs.iter().enumerate() {
                if map[n] == node {
                    return stim.value(i, t);
                }
            }
            0.5 * model.vdd
        }
    };
    let w = transient(&fnet, model, &drive, times, &v0, &[]);
    FlatSim { wave: w, map, net: fnet }
}

/// Energy delivered by the V_dd supply over a transistor-level waveform (J).
pub fn wave_supply_energy(model: &DeviceModel, net: &Net, w: &Wave) -> f64 {
    let mut q = 0.0;
    for k in 1..w.t.len() {
        let i0 = supply_current(net, model, &w.v[k - 1]);
        let i1 = supply_current(net, model, &w.v[k]);
        q += 0.5 * (i0 + i1) * (w.t[k] - w.t[k - 1]);
    }
    model.vdd * q
}
