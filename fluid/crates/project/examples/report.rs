//! Print a project's report: `cargo run -p project --release --example report -- path/to/project.json`.
fn main() {
    let path = std::env::args().nth(1).expect("project path");
    let t = std::time::Instant::now();
    let p = project::Project::load(&path).unwrap_or_else(|e| panic!("{e}"));
    let e = p.evaluate();
    let t1 = t.elapsed();
    print!("{}", project::report::text(&p, &e, true));
    eprintln!("evaluate {:.2?}, with contributions {:.2?}", t1, t.elapsed());
}
