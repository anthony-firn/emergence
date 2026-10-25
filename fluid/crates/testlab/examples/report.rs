//! Print the virtual materials-testing lab report.
//!
//! `cargo run -p testlab --example report`            full standard-test battery
//! `cargo run -p testlab --example report -- --full`  plus MD and plate impact
//!
//! The battery runs every material through the standard tests of the
//! mechanical, thermal, viscoelastic, magnetic, optical, wetting and granular
//! rungs, and `format_report` prints one table per material followed by the
//! cross-material summary of which folds were certified.

fn main() {
    let full = std::env::args().any(|a| a == "--full");
    let rows = if full {
        testlab::report::full_report()
    } else {
        testlab::report::battery_report()
    };
    print!("{}", testlab::report::format_report(&rows));
}
