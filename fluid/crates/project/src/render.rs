//! Pictures of a project's CAD parts, drawn straight from the shared model
//! by the workspace's sphere tracer (`viz`), so the image is the geometry
//! every analysis used, at the project's current values.

use crate::graph::{Evaluation, Project};

/// Render part `id` at its nominal parameter values; PNG bytes.
pub fn part_png(p: &Project, e: &Evaluation, id: &str, width: usize, height: usize) -> Result<Vec<u8>, String> {
    use viz::camera::Camera;
    use viz::scene::{Light, Material, Object, Scene};
    use viz::shade::{render, RenderOptions};

    let part = p.part(id).ok_or_else(|| format!("`{id}` is not a part"))?;
    let run = e.runs.get(id).ok_or_else(|| format!("`{id}` was not evaluated"))?;
    if let Some(err) = &run.error {
        return Err(err.clone());
    }
    let x: Vec<f64> = run.inputs.iter().take(part.params.len()).map(|(_, _, v)| v.map(|v| v.nom).unwrap_or(f64::NAN)).collect();
    let solid = part.solid(&x)?;
    let (lo, hi) = solid.bbox();
    let size = (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
    let centre = [(lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0, (lo[2] + hi[2]) / 2.0];
    // A floor just under the part, checkered at a tenth of its size.
    let floor = cad::Solid::box_between(
        [centre[0] - 4.0 * size, centre[1] - 4.0 * size, lo[2] - 0.02 * size],
        [centre[0] + 4.0 * size, centre[1] + 4.0 * size, lo[2] - 1e-6 * size],
    );
    let mut scene = Scene::new()
        .with_material(Material::Diffuse { albedo: 0.72 })
        .with_material(Material::Checker { a: 0.55, b: 0.42, period: 0.1 * size })
        .with_object(Object::new(id, solid, 0))
        .with_object(Object::new("floor", floor, 1))
        .with_light(Light::Distant { direction: [0.45, 0.65, -0.62], irradiance: 1.0 })
        .with_light(Light::Distant { direction: [-0.7, 0.2, -0.4], irradiance: 0.35 });
    scene.background = 0.06;
    let eye = [centre[0] - 1.7 * size, centre[1] - 2.3 * size, centre[2] + 1.6 * size];
    let camera = Camera::looking_at(eye, centre, 0.42);
    let image = render(&scene, &camera, &RenderOptions::preview(width, height)).beauty;
    let tmp = std::env::temp_dir().join(format!("fluid-project-render-{}-{id}.png", std::process::id()));
    image.write_png(&tmp, 0.26).map_err(|e| e.to_string())?;
    let bytes = std::fs::read(&tmp).map_err(|e| e.to_string())?;
    let _ = std::fs::remove_file(&tmp);
    Ok(bytes)
}

/// Standard base64, for embedding a PNG in a self-contained page.
pub fn base64(bytes: &[u8]) -> String {
    const T: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut s = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for c in bytes.chunks(3) {
        let n = (c[0] as u32) << 16 | (*c.get(1).unwrap_or(&0) as u32) << 8 | *c.get(2).unwrap_or(&0) as u32;
        for k in 0..4 {
            if k <= c.len() {
                s.push(T[(n >> (18 - 6 * k) & 63) as usize] as char);
            } else {
                s.push('=');
            }
        }
    }
    s
}
