//! Conservative SI adapter for explicitly scheduled kinetic-particle sources.
use crate::{Solver, System};
use standard_model::{coupling, FourMomentum};

/// Apply (interior cell, transferred total four-momentum in GeV) sources.
/// The interaction closure must determine what transfers to matter: neutrinos
/// and unresolved hadrons are not automatically absorbed. The existing fluid
/// source checks every cell before committing and includes momentum work.
pub fn deposit<const N: usize, S: System<N>>(
    solver: &mut Solver<N, S>,
    transfers: &[(usize, FourMomentum)],
) -> Result<f64, &'static str> {
    let tallies: Result<Vec<_>, _> = transfers
        .iter()
        .map(|(cell, transfer)| coupling::deposit(*cell, *transfer))
        .collect();
    crate::radiation::deposit(solver, &tallies?)
}
