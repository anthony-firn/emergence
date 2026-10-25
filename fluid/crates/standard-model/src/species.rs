//! PDG 2025 central masses/widths. Quark masses are running mass parameters,
//! not measured free-particle pole masses (top uses the quoted direct mass).
use crate::HBAR_GEV_S;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Domain {
    /// Colored final states are retained as unresolved hadronic events.
    Vacuum,
    /// Color-averaged on-shell partons: no confinement, shower, or QCD forces.
    Partonic,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Species {
    Down,
    Up,
    Strange,
    Charm,
    Bottom,
    Top,
    AntiDown,
    AntiUp,
    AntiStrange,
    AntiCharm,
    AntiBottom,
    AntiTop,
    Electron,
    Positron,
    Muon,
    AntiMuon,
    Tau,
    AntiTau,
    ElectronNeutrino,
    AntiElectronNeutrino,
    MuonNeutrino,
    AntiMuonNeutrino,
    TauNeutrino,
    AntiTauNeutrino,
    Photon,
    Gluon0,
    Gluon1,
    Gluon2,
    Gluon3,
    Gluon4,
    Gluon5,
    Gluon6,
    Gluon7,
    WPlus,
    WMinus,
    Z,
    Higgs,
}
use Species::*;
pub const ALL: [Species; 37] = [
    Down,
    Up,
    Strange,
    Charm,
    Bottom,
    Top,
    AntiDown,
    AntiUp,
    AntiStrange,
    AntiCharm,
    AntiBottom,
    AntiTop,
    Electron,
    Positron,
    Muon,
    AntiMuon,
    Tau,
    AntiTau,
    ElectronNeutrino,
    AntiElectronNeutrino,
    MuonNeutrino,
    AntiMuonNeutrino,
    TauNeutrino,
    AntiTauNeutrino,
    Photon,
    Gluon0,
    Gluon1,
    Gluon2,
    Gluon3,
    Gluon4,
    Gluon5,
    Gluon6,
    Gluon7,
    WPlus,
    WMinus,
    Z,
    Higgs,
];
pub const GLUONS: [Species; 8] = [
    Gluon0, Gluon1, Gluon2, Gluon3, Gluon4, Gluon5, Gluon6, Gluon7,
];
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Charges {
    /// Electric charge in thirds of e, baryon number in thirds of B.
    pub electric3: i32,
    pub baryon3: i32,
    /// Electron, muon, tau lepton numbers (massless, no mixing closure).
    pub lepton: [i32; 3],
}
impl std::ops::Add for Charges {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            electric3: self.electric3 + b.electric3,
            baryon3: self.baryon3 + b.baryon3,
            lepton: std::array::from_fn(|i| self.lepton[i] + b.lepton[i]),
        }
    }
}
impl Species {
    pub fn conjugate(self) -> Self {
        match self {
            Down => AntiDown,
            Up => AntiUp,
            Strange => AntiStrange,
            Charm => AntiCharm,
            Bottom => AntiBottom,
            Top => AntiTop,
            AntiDown => Down,
            AntiUp => Up,
            AntiStrange => Strange,
            AntiCharm => Charm,
            AntiBottom => Bottom,
            AntiTop => Top,
            Electron => Positron,
            Positron => Electron,
            Muon => AntiMuon,
            AntiMuon => Muon,
            Tau => AntiTau,
            AntiTau => Tau,
            ElectronNeutrino => AntiElectronNeutrino,
            AntiElectronNeutrino => ElectronNeutrino,
            MuonNeutrino => AntiMuonNeutrino,
            AntiMuonNeutrino => MuonNeutrino,
            TauNeutrino => AntiTauNeutrino,
            AntiTauNeutrino => TauNeutrino,
            WPlus => WMinus,
            WMinus => WPlus,
            _ => self,
        }
    }
    pub fn pdg_id(self) -> i32 {
        match self {
            Down => 1,
            Up => 2,
            Strange => 3,
            Charm => 4,
            Bottom => 5,
            Top => 6,
            AntiDown => -1,
            AntiUp => -2,
            AntiStrange => -3,
            AntiCharm => -4,
            AntiBottom => -5,
            AntiTop => -6,
            Electron => 11,
            Positron => -11,
            Muon => 13,
            AntiMuon => -13,
            Tau => 15,
            AntiTau => -15,
            ElectronNeutrino => 12,
            AntiElectronNeutrino => -12,
            MuonNeutrino => 14,
            AntiMuonNeutrino => -14,
            TauNeutrino => 16,
            AntiTauNeutrino => -16,
            Photon => 22,
            WPlus => 24,
            WMinus => -24,
            Z => 23,
            Higgs => 25,
            _ => 21,
        }
    }
    /// PDG 21 is color inclusive; choose Gluon0, or select a GLUONS entry.
    pub fn from_pdg(id: i32) -> Option<Self> {
        ALL.into_iter().find(|s| s.pdg_id() == id)
    }
    pub fn fermion(self) -> bool {
        self.pdg_id().abs() <= 16
    }
    pub fn twice_spin(self) -> u8 {
        if self.fermion() {
            1
        } else if self == Higgs {
            0
        } else {
            2
        }
    }
    pub fn quark(self) -> bool {
        self.pdg_id().abs() <= 6
    }
    pub fn colored(self) -> bool {
        self.quark() || self.pdg_id() == 21
    }
    pub fn neutrino(self) -> bool {
        matches!(self.pdg_id().abs(), 12 | 14 | 16)
    }
    pub fn mass(self) -> f64 {
        match self.pdg_id().abs() {
            1 => 0.00470,
            2 => 0.00216,
            3 => 0.0935,
            4 => 1.273,
            5 => 4.183,
            6 => 172.56,
            11 => 0.000_510_998_950,
            13 => 0.105_658_375_5,
            15 => 1.77693,
            24 => 80.369,
            23 => 91.1880,
            25 => 125.20,
            // Minimal SM massless flavor transport; no guessed neutrino masses.
            _ => 0.0,
        }
    }
    pub fn width(self) -> f64 {
        match self.pdg_id().abs() {
            6 => 1.42,
            13 => 2.995_983_6e-19,
            15 => 2.267e-12,
            24 => 2.14,
            23 => 2.4955,
            25 => 3.7e-3,
            _ => 0.0,
        }
    }
    /// None for light quarks is an absent free-parton decay closure, not stability.
    pub fn lifetime(self) -> Option<f64> {
        let w = self.width();
        if w > 0.0 {
            Some(HBAR_GEV_S / w)
        } else {
            None
        }
    }
    pub fn needs_qcd_closure(self) -> bool {
        self.colored() && self.pdg_id().abs() != 6
    }
    pub fn charges(self) -> Charges {
        let id = self.pdg_id();
        let sign = if id < 0 { -1 } else { 1 };
        let mut q = Charges::default();
        if self.quark() {
            q.baryon3 = sign;
            q.electric3 = sign * if id.abs() % 2 == 0 { 2 } else { -1 };
        } else if self.fermion() {
            q.lepton[((id.abs() - 11) / 2) as usize] = sign;
            if !self.neutrino() {
                q.electric3 = -3 * sign;
            }
        } else if id.abs() == 24 {
            q.electric3 = 3 * sign;
        }
        q
    }
    /// Includes 3 quark colors, but excludes antiparticles (explicit entries).
    /// Each gluon entry counts 2 helicities; massless neutrinos count 1.
    pub fn degeneracy(self) -> f64 {
        if self.quark() {
            6.0
        } else if self.neutrino() {
            1.0
        } else if self.fermion() || self.mass() == 0.0 {
            2.0
        } else if self == Higgs {
            1.0
        } else {
            3.0
        }
    }
    pub fn allowed(self, domain: Domain) -> bool {
        domain == Domain::Partonic || !self.needs_qcd_closure()
    }
}
