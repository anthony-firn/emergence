//! Time one standard-cell library build (`cargo run -p semiconductor --release --example libtime`).
use semiconductor::gates::{DeviceModel, Ekv};
use semiconductor::timing::build_library;
fn main() {
    let (n, p) = (Ekv::n(), Ekv::p());
    let m = DeviceModel::new(&n, &p);
    let t = std::time::Instant::now();
    let lib = build_library(&m, &[15e-12, 40e-12, 120e-12, 350e-12], &[1e-15, 4e-15, 12e-15], &[(55e-12, 2e-15), (55e-12, 8e-15), (180e-12, 2e-15), (180e-12, 8e-15)]);
    let fa = &lib.cells[4];
    println!("library in {:.2} s; FA sum delay {:.6e}, cert {:.4}", t.elapsed().as_secs_f64(), fa.outs[0][0].rise.interp(100e-12, 4e-15).0, fa.outs[0][0].rise.cert_error);
}
