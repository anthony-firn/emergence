//! Fixed source/query set; direct and tree both validate sources once.
use gravity::{norm, Body, Direct, Tree, Vec3, G};
use std::{hint::black_box, time::Instant};

fn main() -> Result<(), &'static str> {
    let n = 32768;
    let queries = 256;
    let radius = 6.957e8;
    let mut state = 942_u64;
    let mut uniform = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 * (1.0 / (1_u64 << 53) as f64)
    };
    let bodies: Vec<_> = (0..n)
        .map(|_| Body {
            position: Vec3(std::array::from_fn(|_| radius * (2.0 * uniform() - 1.0))),
            mass: 1.3271244e20 / (G * n as f64),
        })
        .collect();
    let epsilon = 0.01 * radius;
    let theta = 0.35;
    let start = Instant::now();
    let tree = Tree::new(bodies.clone(), epsilon, theta)?;
    let build = start.elapsed().as_secs_f64();
    let direct = Direct::new(bodies.clone(), epsilon)?;
    // Warm both paths; timings exclude construction and source validation.
    black_box(direct.evaluate(bodies[0].position, Some(0))?);
    black_box(tree.evaluate(bodies[0].position, Some(0))?);
    let start = Instant::now();
    let reference: Vec<_> = (0..queries)
        .map(|k| {
            let i = k * n / queries;
            direct.evaluate(bodies[i].position, Some(i))
        })
        .collect::<Result<_, _>>()?;
    let direct_seconds = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let actual: Vec<_> = (0..queries)
        .map(|k| {
            let i = k * n / queries;
            tree.evaluate(bodies[i].position, Some(i))
        })
        .collect::<Result<_, _>>()?;
    let tree_seconds = start.elapsed().as_secs_f64();
    let error: f64 = actual
        .iter()
        .zip(&reference)
        .map(|((a, _), b)| norm(a.acceleration - b.acceleration).powi(2))
        .sum();
    let power: f64 = reference.iter().map(|b| norm(b.acceleration).powi(2)).sum();
    let pairs: usize = actual.iter().map(|(_, work)| work.direct_pairs).sum();
    println!("synthetic uniform cube: sources={n} queries={queries} theta={theta} epsilon/R=0.01");
    println!(
        "tree_build_s={build:.6} direct_query_s={direct_seconds:.6} tree_query_s={tree_seconds:.6}"
    );
    println!(
        "query_speedup={:.3} acceleration_relative_rms={:.9e}",
        direct_seconds / tree_seconds,
        (error / power).sqrt()
    );
    println!("direct_pairs={pairs} reference_pairs={}", queries * (n - 1));
    Ok(())
}
