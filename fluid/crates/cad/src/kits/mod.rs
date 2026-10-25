//! Generators: the parametric things a shop, a fab or a mill can make.
//!
//! Each kit turns a handful of engineering parameters into a [`Solid`] plus
//! the non-geometric data that goes with it — a layer stack, a strut graph,
//! a doping profile, a yarn count. Each one is validated against an exact
//! identity or a published standard; see `analysis/cad.md`.
//!
//! * [`part`] — machined and printed solids: plates, bosses, ribs, shafts,
//!   involute spur gears, ISO metric threads.
//! * [`lattice`] — periodic strut lattices and triply periodic minimal
//!   surfaces, with the strut graph `architected` homogenises.
//! * [`pcb`] — printed circuit boards: stackup, copper polygons, traces,
//!   pads, vias, and the impedance and resistance of what was drawn.
//! * [`silicon`] — a process flow (substrate, deposit, etch, implant) whose
//!   output is the region map and doping closure `semiconductor` takes.
//! * [`textile`] — woven and weft-knitted fabrics: yarn paths, swept yarn
//!   solids, cover factor and fibre volume fraction.
//! * [`pipe`] — pipe and duct runs: standard schedules, elbows, reducers
//!   and tees, the bore as a fluid domain, and the hydraulics of what was
//!   drawn.
//! * [`section`] — structural shapes and the section properties a beam
//!   element needs.
//! * [`airfoil`] — NACA four-digit sections and lofted blades, the one
//!   family here defined by a published polynomial rather than a drawing.
//! * [`optic`] — lenses, mirrors, prisms and windows from their sag
//!   equations, which `waveoptics` then traces real rays through.
//! * [`crystal`] — atomistic structures: crystals, nanoparticles cut from
//!   them by any solid, graphene sheets and rolled nanotubes.
//! * [`sheet`] — formed sheet-metal parts and the flat patterns they are
//!   cut from, with the standard bend-allowance bookkeeping.

use serde_json::Value as Json;

use crate::solid::Solid;

pub mod airfoil;
pub mod crystal;
pub mod lattice;
pub mod part;
pub mod optic;
pub mod pcb;
pub mod pipe;
pub mod section;
pub mod sheet;
pub mod silicon;
pub mod textile;

/// Build a kit solid from `{"op":"kit","kit":"<name>","args":{...}}` in a
/// JSON model. The argument structs are the kits' own `*Spec` types.
pub fn from_json(kit: &str, args: &Json) -> Result<Solid, String> {
    fn parse<T: serde::de::DeserializeOwned>(kit: &str, args: &Json) -> Result<T, String> {
        serde_json::from_value(args.clone()).map_err(|e| format!("kit {kit}: {e}"))
    }
    match kit {
        "spur_gear" => Ok(parse::<part::SpurGear>(kit, args)?.solid()),
        "thread" => Ok(parse::<part::Thread>(kit, args)?.solid()),
        "plate" => Ok(parse::<part::Plate>(kit, args)?.solid()),
        "strut_lattice" => {
            let spec: lattice::LatticeSpec = parse(kit, args)?;
            Ok(spec.graph().solid(spec.counts))
        }
        "tpms" => Ok(parse::<lattice::TpmsSpec>(kit, args)?.solid()),
        "woven" => Ok(parse::<textile::Weave>(kit, args)?.solid()),
        "knit" => Ok(parse::<textile::Knit>(kit, args)?.solid()),
        "spring" => Ok(parse::<part::Spring>(kit, args)?.solid()?.0),
        "beam" => {
            #[derive(serde::Deserialize)]
            struct Spec {
                section: section::Section,
                length: f64,
            }
            let spec: Spec = parse(kit, args)?;
            spec.section.beam(spec.length)
        }
        "wing" => {
            #[derive(serde::Deserialize)]
            struct Spec {
                #[serde(flatten)]
                section: airfoil::Naca4,
                span: f64,
            }
            let spec: Spec = parse(kit, args)?;
            Ok(spec.section.wing(spec.span))
        }
        "lens" => parse::<optic::Lens>(kit, args)?.solid(),
        "mirror" => parse::<optic::Mirror>(kit, args)?.solid(),
        "prism" => Ok(parse::<optic::Prism>(kit, args)?.solid()),
        "window" => Ok(parse::<optic::Window>(kit, args)?.solid()),
        "sheet_part" => {
            #[derive(serde::Deserialize)]
            struct Spec {
                #[serde(flatten)]
                part: sheet::SheetPart,
                #[serde(default = "sheet_points")]
                per_bend: usize,
            }
            fn sheet_points() -> usize {
                24
            }
            let spec: Spec = parse(kit, args)?;
            spec.part.solid(spec.per_bend)
        }
        "pipe_run" => {
            #[derive(serde::Deserialize)]
            struct Spec {
                pipe: pipe::Pipe,
                run: pipe::Run,
                /// `bore` for the fluid domain, `wall` for the metal.
                #[serde(default = "pipe_which")]
                part: String,
                #[serde(default = "pipe_sides")]
                sides: usize,
            }
            fn pipe_which() -> String {
                "wall".into()
            }
            fn pipe_sides() -> usize {
                48
            }
            let spec: Spec = parse(kit, args)?;
            let path = spec.run.centreline()?;
            Ok(match spec.part.as_str() {
                "bore" => spec.pipe.bore(&path, spec.sides).0,
                "wall" => spec.pipe.wall_solid(&path, spec.sides).0,
                other => return Err(format!("pipe_run: part must be \"bore\" or \"wall\", got {other:?}")),
            })
        }
        other => Err(format!(
            "unknown kit {other:?} (spur_gear, thread, plate, strut_lattice, tpms, woven, knit, \
             spring, beam, wing, lens, mirror, prism, window, sheet_part, pipe_run)"
        )),
    }
}
