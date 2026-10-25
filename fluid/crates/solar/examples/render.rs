//! Render the actual GONG-initialized surface and its computed PFSS topology.
mod common;
use gravity::{norm, Vec3};
use rayon::prelude::*;
use solar::{
    corona::Pfss,
    photosphere,
    render::{self, Options as RenderOptions, View},
    surface::Parameters,
    trace::{self, Options},
};
use std::{
    f64::consts::TAU,
    fs,
    io::{BufWriter, Write},
    path::PathBuf,
    time::Instant,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let output = PathBuf::from(
        arguments
            .first()
            .map(String::as_str)
            .unwrap_or("/tmp/solar-render"),
    );
    let frames: usize = arguments
        .get(1)
        .map(|v| v.parse())
        .transpose()?
        .unwrap_or(28);
    let orbit = arguments.get(2).is_some_and(|s| s == "orbit");
    if frames == 0 || frames > 1000 {
        return Err("frame count must be 1..=1000".into());
    }
    fs::create_dir_all(&output)?;
    let horizon_days = if frames == 1 { 0.0 } else { 27.0 };
    let mut surface = common::map(common::INITIAL, Parameters::default())?;
    let removed_mean = surface.balance_flux();
    let initial_flux = surface.net_flux();
    let teff = (solar::L_SUN
        / (4.0 * std::f64::consts::PI * solar::R_SUN.powi(2) * photosphere::SIGMA))
        .sqrt()
        .sqrt();
    let mut options = RenderOptions {
        width: 900,
        height: 900,
        samples_per_axis: 2,
        view: View {
            longitude: 220_f64.to_radians(),
            latitude: 7_f64.to_radians(),
            extent: 1.12,
        },
    };
    let start = Instant::now();
    let continuum = render::continuum(teff, options)?;
    continuum.write_pfm(output.join("continuum.pfm"))?;
    let peak = photosphere::rgb(1.0, teff)?.into_iter().fold(0.0, f64::max);
    continuum.write_png(output.join("continuum.png"), 2.0 / peak)?;
    // Fixed equal-area geometry, independent of later maps, avoids flickering
    // random seeds. Curves are not weighted as coronal emissivity or open flux.
    let seeds: Vec<_> = (0..10)
        .flat_map(|i| {
            (0..24).map(move |j| {
                let mu = -1.0 + (i as f64 + 0.5) * 0.2;
                let phi = (j as f64 + 0.5) * TAU / 24.0;
                let sin = (1.0 - mu * mu).sqrt();
                Vec3([sin * phi.cos(), sin * phi.sin(), mu]) * (solar::R_SUN * 1.0001)
            })
        })
        .collect();
    let mut manifest = BufWriter::new(fs::File::create(output.join("frames.csv"))?);
    writeln!(manifest,"frame,days,view_longitude_deg,closed,open,unresolved,dipole_G,signed_flux_drift_Wb,open_flux_Wb,runtime_s")?;
    for frame in 0..frames {
        let frame_start = Instant::now();
        let days = if frames == 1 {
            0.0
        } else {
            27.0 * frame as f64 / (frames - 1) as f64
        };
        surface.advance(days * 86400.0 - surface.time())?;
        let field = Pfss::from_surface(&surface, 10, 2.5)?;
        let lines = seeds
            .par_iter()
            .map(|p| trace::trace(&field, *p, Options::default()))
            .collect::<Result<Vec<_>, _>>()?;
        let closed = lines.iter().filter(|l| l.closed()).count();
        let open = lines.iter().filter(|l| l.open()).count();
        let unresolved = lines.len() - closed - open;
        let longitude = 220.0
            - if orbit && frames > 1 {
                360.0 * frame as f64 / (frames - 1) as f64
            } else {
                0.0
            };
        options.view = View {
            longitude: longitude.to_radians(),
            latitude: 7_f64.to_radians(),
            extent: 2.65,
        };
        render::magnetic(&surface, &lines, 20e-4, options)?
            .write_png(output.join(format!("frame{frame:04}.png")), 2.0)?;
        // Raw SI geometry and explicit boundary outcomes, not a drawn-loop fit.
        let mut geometry = BufWriter::new(fs::File::create(
            output.join(format!("trace{frame:04}.csv")),
        )?);
        writeln!(geometry, "line,vertex,x_m,y_m,z_m,end_a,end_b")?;
        for (index, line) in lines.iter().enumerate() {
            for (vertex, p) in line.points.iter().enumerate() {
                writeln!(
                    geometry,
                    "{index},{vertex},{:.12e},{:.12e},{:.12e},{},{}",
                    p.0[0],
                    p.0[1],
                    p.0[2],
                    line.ends[0].label(),
                    line.ends[1].label()
                )?;
                if norm(*p) < solar::R_SUN * (1.0 - 1e-12) {
                    return Err("trace entered opaque photosphere".into());
                }
            }
        }
        let mut map = BufWriter::new(fs::File::create(output.join(format!("map{frame:04}.csv")))?);
        for row in surface.field().chunks(surface.dimensions().1) {
            writeln!(
                map,
                "{}",
                row.iter()
                    .map(|v| format!("{v:.12e}"))
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
        }
        let open_flux = field.open_flux(60, 120)?;
        let elapsed = frame_start.elapsed().as_secs_f64();
        writeln!(manifest,"{frame},{days:.9},{longitude:.9},{closed},{open},{unresolved},{:.9},{:.9e},{open_flux:.9e},{elapsed:.6}",surface.axial_dipole()*1e4,surface.net_flux()-initial_flux)?;
        println!("frame={frame}/{frames} day={days:.3} closed={closed} open={open} unresolved={unresolved} elapsed_s={elapsed:.3}");
    }
    let metadata = format!(
        r#"{{
  "initial_observation": "GONG 2019-01-08T11:14 UTC",
  "model": "prescribed surface flux transport, no new emergence; frozen PFSS per frame",
  "diffusion": "L-stable SDIRK2, cached tridiagonal solves",
  "continuum": "Eddington grey plane-parallel LTE, formal spectral transfer, CIE visible response",
  "effective_temperature_K": {teff:.9},
  "continuum_png_exposure": {exposure:.12e},
  "continuum_pfm": "linear RGB continuum radiance; no display exposure",
  "magnetic_colors": "diagnostic asinh scale +/-20 G; gold=closed, green=open, grey=unresolved; no emission",
  "camera": "{camera}; fixed latitude 7 degrees; straight orthographic rays, no ephemeris",
  "seed_count": {seed_count},
  "trace_tolerance_R": 1e-6,
  "source_surface_R": 2.5,
  "harmonic_degree": 10,
  "initial_removed_mean_G": {removed:.9},
  "frames": {frames},
  "days": {horizon_days},
  "runtime_s": {runtime:.6}
}}
"#,
        exposure = 2.0 / peak,
        camera = if orbit {
            "inspection turntable, not solar rotation"
        } else {
            "fixed Carrington longitude 220 degrees"
        },
        seed_count = seeds.len(),
        removed = removed_mean * 1e4,
        runtime = start.elapsed().as_secs_f64()
    );
    fs::write(output.join("metadata.json"), metadata)?;
    let html = include_str!("render_viewer.html")
        .replace("__FRAMES__", &frames.to_string())
        .replace(
            "__CAMERA__",
            if orbit {
                "Inspection camera rotates; this is not solar rotation."
            } else {
                "Fixed Carrington camera."
            },
        );
    fs::write(output.join("index.html"), html)?;
    println!(
        "output={} runtime_s={:.3} Teff_K={teff:.3}",
        output.display(),
        start.elapsed().as_secs_f64()
    );
    Ok(())
}
