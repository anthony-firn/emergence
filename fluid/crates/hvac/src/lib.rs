//! Heating, ventilation, air conditioning and cooling design.
//!
//! Each module is the standard engineering method, implemented as
//! published, and checked against an exact or published result:
//! psychrometrics (ASHRAE Fundamentals Ch. 1), thermal comfort (ISO 7730).

pub mod coil;
pub mod comfort;
pub mod envelope;
pub mod hydraulics;
pub mod psychro;
pub mod refrigeration;
pub mod solar;
pub mod zone;
