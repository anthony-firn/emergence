//! Structural reference, exterior GR photon and efficient polar-field transport.
mod common;
use gravity::{
    metric::{Metric, Schwarzschild},
    norm, Field, Vec3,
};
use solar::{corona::Pfss, structure::Profile, surface::Parameters, GM_SUN, R_SUN};
use standard_model::{
    geodesic::{propagate, Options},
    FourMomentum,
};
use std::time::Instant;

fn main() -> Result<(), &'static str> {
    let profile = Profile::model_s()?;
    let field = profile.gravity(Vec3::default())?;
    println!(
        "published Model S: knots={} mass_kg={:.9e} outer_table_radius_m={:.9e}",
        profile.layers.len(),
        field.total_mass(),
        field.radius()
    );
    println!(
        "center_density_kg_m3={:.9e} center_temperature_K={:.9e} surface_gravity_m_s2={:.9e}",
        profile.layers[0].density,
        profile.layers[0].temperature,
        norm(field.sample(Vec3([field.radius(), 0.0, 0.0]))?.acceleration)
    );

    // Vacuum outside a nominal spherical Sun, separate from the tabulated star.
    let metric = Schwarzschild {
        gm: GM_SUN,
        center: Vec3::default(),
    };
    let x = Vec3([R_SUN, 0.0, 0.0]);
    let photon = FourMomentum::on_shell(0.0, Vec3([1e-9, 0.0, 0.0]))?;
    let orbit = propagate(
        &metric,
        x,
        photon,
        0.0,
        100.0,
        f64::INFINITY,
        Options::default(),
    )?;
    let initial = metric.geometry(x)?.lapse * photon.energy;
    let final_energy = metric.geometry(orbit.position)?.lapse * orbit.four.energy;
    println!(
        "exterior_photon_100s: local_energy_ratio={:.12} killing_energy_relative_error={:.9e}",
        orbit.four.energy / photon.energy,
        (final_energy - initial) / initial
    );

    let mut surface = common::map(common::INITIAL, Parameters::default())?;
    surface.balance_flux();
    let corona = Pfss::from_surface(&surface, 8, 2.5)?;
    println!(
        "initial_PFSS: lmax=8 source_surface/R=2.5 open_flux_Wb={:.9e}",
        corona.open_flux(120, 240)?
    );
    let initial_flux = surface.net_flux();
    let mut folded = surface.axisymmetric()?;
    let start = Instant::now();
    let steps = surface.advance(365.25 * 86400.0)?;
    let full_seconds = start.elapsed().as_secs_f64();
    let start = Instant::now();
    folded.advance(365.25 * 86400.0)?;
    let folded_seconds = start.elapsed().as_secs_f64();
    let max_mean_error = surface
        .mean_field()
        .iter()
        .zip(folded.field())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!("conditional_1yr: grid=60x120 steps={steps} full_s={full_seconds:.6} folded_60x1_s={folded_seconds:.6} speedup={:.3}", full_seconds / folded_seconds);
    println!("fold_max_mean_difference_T={max_mean_error:.9e} full_dipole_G={:.9} folded_dipole_G={:.9} signed_flux_drift_Wb={:.9e}", surface.axial_dipole()*1e4, folded.axial_dipole()*1e4, surface.net_flux()-initial_flux);
    println!("No new emergences: this is conditional transport, not a predicted solar cycle or eruption.");
    Ok(())
}
