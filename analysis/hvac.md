# HVAC: heating, cooling and ventilation design

Crate `fluid/crates/hvac`; refrigerants in `fluid/crates/compressible`;
project analyses in `fluid/crates/project` (`zone_loads`,
`refrigeration_cycle`, `duct_run`, `fan`, `cooling_coil`,
`thermal_comfort`); example `fluid/projects/office-hvac/`.

Each module is a standard engineering method implemented as published, and
each is checked against an exact result or an independent reference. Where
a published value could not be confirmed, or a published accuracy claim did
not hold, this document says so.

## Refrigerants

R-32, R-134a, R-410A (pseudo-pure), R-290, R-1234yf and R-744 (CO₂) run on
the repository's Helmholtz equation-of-state machinery. The coefficients
come from the fluid files bundled with CoolProp (`data/build_refrigerants.py`,
read locally), which also records CoolProp's own values for the tests.
CoolProp's reference-state offset a₁ + a₂τ maps onto two existing terms.

| check against CoolProp | result |
|---|---|
| pressure, single phase | ≤ 1.1e-15 |
| enthalpy, entropy | ≤ 5.8e-10 |
| c_p, sound speed (second derivatives) | ≤ 1.4e-9 |
| saturation pressure and densities | ≤ 5e-11 |
| R-410A: the one Maxwell pressure lies between CoolProp's bubble and dew pressures | at every temperature (e.g. 1429.9 < 1431.1 < 1434.3 kPa at 292.8 K) |

Tolerances are set by derivative order (1e-9 for first derivatives, 1e-8
for second). A first blanket 1e-9 was too tight for R-32's c_p (1.4e-9).

## Moist air and comfort

**Psychrometrics** follow ASHRAE Fundamentals 2017 Ch. 1:

* Hyland–Wexler saturation over water and over ice;
* humidity ratio, enthalpy and specific volume;
* dew point, and wet bulb from the psychrometric equation;
* pressure at altitude.

| check | result |
|---|---|
| saturation vs IAPWS-95 (0.01–100 °C) | within 2.2e-4 |
| wet bulb = dew point = dry bulb at saturation | to 1e-6 K |
| ASHRAE Example 2 (40 °C / 20 °C wet bulb → W = 0.0065) | reproduced (0.00646) with the example's tabulated W_s, which includes the enhancement factor; the equations without it read 0.9 % low there, the stated error of omitting it |
| standard pressure at 1500 m | 84.556 kPa to 2e-4 |

**Comfort** is Fanger's PMV/PPD, computed by ISO 7730 Annex D's own
program, plus operative temperature. Ten of the thirteen rows of Table D.1
(as recalled) reproduce to the table's two decimals. Three recalled rows do
not, and could not be checked against the standard's text offline. They are
printed, not asserted. Monotonicity is checked: PMV rises with air
temperature, radiant temperature, activity and clothing, and falls with air
speed.

## Envelope

**Steady:** U from layers, with ISO 6946 surface resistances.

**Periodic:** ISO 13786 heat-transfer matrices give the decrement factor,
time lag and internal admittance. The per-layer matrix was checked term by
term against the standard's closed forms. A first version multiplied the
layers in reverse order (which leaves Y₁₂ unchanged but swaps the inside and
outside admittances) and carried a stray −π in the lag. Both were fixed to
the standard's ordering.

**Transient:** a finite-volume wall, BDF2 in time. Crank–Nicolson was tried
first: at hourly steps it never damps the stiff modes of thin layers, so
steady state was never reached.

| check | result |
|---|---|
| finite-volume wall vs exact periodic response, 2/4/8 nodes per layer | second order (error ratios 3.71, 3.93); Richardson within 9e-5; lag 9.819 vs 9.818 h |
| steady state | U·ΔT to 1e-9 |
| very slow swings | decrement → 1; the lag tends to the wall's finite delay, not to zero |
| masonry wall, 80 mm insulation | U 0.365 W/m²K, decrement 0.263, lag 9.8 h, admittance 4.9 W/m²K |

## Sun and loads

**Solar** follows ASHRAE Fundamentals Ch. 14 as published:

* equation of time, declination, solar time;
* altitude and azimuth (atan2: arccos lost 1e-6° at noon);
* the clear-sky model with site optical depths and Kasten–Young air mass;
* beam, diffuse (ASHRAE's Y factor on verticals) and ground-reflected
  irradiance on any surface;
* sol-air temperature and the design-day temperature profile.

**Zone loads** come from a heat balance, with the room held at its set
point:

* transient walls driven by sol-air temperature;
* internal elements as half-elements with an adiabatic mid-plane;
* windows by U and SHGC;
* radiant gains absorbed by the surfaces by area, convective gains at once;
* outdoor air sensible and latent.

Cooling uses a periodic design day; heating is steady at the winter design
temperature.

| check | result |
|---|---|
| solar noon altitude 90° − L + δ; sunrise at cos H = −tan L tan δ | to 1e-9 |
| ASHRAE's declination (Cooper) vs Spencer's series | within 1.39° all year (the formula's accuracy) |
| no sun, constant outdoor air, no gains | load = UA·ΔT + ventilation, exactly |
| periodic design day: mean load = steady balance of the mean drivers | to 3.4e-9 (exact for the linear BDF2 system; includes the 1 − U·R_si share of absorbed radiant heat that returns to the room) |
| 1 kW of gain 8–18 h, radiant vs convective | radiant peaks at 791 W at 17:45; convective at 1000 W at 08:00 |

A note on test design: the first "no sun" test used a huge optical depth.
That does not darken the ASHRAE sky, because the diffuse exponent turns
negative. The test now puts the room in polar night.

## Distribution and plant

**Ducts and pipes:** Darcy–Weisbach. The laminar factor is exact; Colebrook
is solved by Newton. The transitional band (Re 2300–4000) is flagged and
given the turbulent value. Rectangular ducts use ASHRAE's equivalent
diameter.

**Fans and pumps:** a quadratic through three catalogue points, and the
operating point where it meets the system curve.

| check | result |
|---|---|
| laminar drop vs Hagen–Poiseuille | 1e-12 |
| Colebrook residual | 1e-12 |
| Haaland vs Colebrook (published ~1.5 %) | 1.34 % |
| Swamee–Jain vs Colebrook (often quoted at 1 %) | **2.83 %** at its domain's corner (Re 5000, ε/D 0.01) |
| operating point on a quadratic system vs the exact root | 1e-12; a system the fan cannot drive is reported, not extrapolated |

**Vapour-compression cycle:** superheat, isentropic efficiency, subcooling,
isenthalpic expansion. Checked against the same cycle computed by CoolProp
(5/45 °C, 5 K superheat, 3 K subcooling, η 0.7):

| refrigerant | COP | vs CoolProp |
|---|---|---|
| R-134a | 4.00765 | 2.4e-9 |
| R-32 | 3.80451 | 1.1e-9 |
| R-290 | 3.93318 | 1.2e-10 |
| R-1234yf | 3.84115 | 7e-12 |
| R-410A | 3.692 | +0.28 % (its glide: CoolProp uses separate dew and bubble pressures) |

The COP tolerance is derived: the enthalpies agree to 1e-9, and the COP
divides by the compressor work, which amplifies that by h₂/(h₂ − h₁).
Every COP is below Carnot. Transcritical CO₂ is refused with the reason.

**Coils:** apparatus dew point and bypass factor, dry when appropriate, and
ε-NTU for sensible exchangers. The total equals the enthalpy drop, and the
latent part equals the condensate times its latent heat, to 1e-9. Every
arrangement reduces to 1 − e^−NTU at C_r = 0, and crossflow lies between
parallel and counterflow.

## In the project graph: the corner office

`projects/office-hvac` is the south-west corner office (120 m²) of a New
York floor, owned by five teams:

* architecture: facade;
* electrical: lighting, plug loads, panel;
* mechanical: ventilation, duct, fan, R-32 plant;
* systems: design climate, marked as estimates to check against the site's
  ASHRAE data;
* the owner: comfort and efficiency.

| result | value |
|---|---|
| peak cooling | 7.96 kW at 15:30 (6.25 sensible + 1.71 latent; 66 W/m²), 4.0 kW of it solar through the glass |
| heating | 4.04 kW |
| supply air | 469 L/s (3.9 L/s/m²); duct 3.5 m/s |
| plant | R-32 COP 4.07, compressor 1.96 kW |
| comfort | PMV −0.24 |
| heating capacity 4 kW ≥ load | **undecided**: outdoor-air tolerance 72 %, glazing U 28 % |
| fan delivers the design airflow | **undecided**: AHU pressure-drop estimate, plug loads, glazing SHGC |

A change to the glazing SHGC moves the solar gain, the cooling load, the
plant's capacity and compressor power, the fan's operating point and the
electrical demand. The report names which team's input each uncertainty
comes from. The building evaluates in 26 ms, 90 ms with contributions, so it
is interactive in `fluid project serve`.

## Honest limits

* **Design-day methods, not annual energy.** There is no 8760-hour
  simulation with weather files, no part-load equipment curves, and no
  energy cost.
* **Room-level heat balance with combined surface coefficients.** There is
  no surface-to-surface long-wave exchange, no shading, and no angular
  SHGC. It is not yet compared with BESTEST (ASHRAE Standard 140).
* **No room airflow.** The CFD here is laminar, so draughts, stratification
  and diffuser throw are out of reach. Comfort assumes a uniform room.
* **Equipment.** Single-stage cycles at design conditions, with no
  heat-exchanger sizing for the refrigerant side. Coils use the ADP/bypass
  model, not a detailed fin-tube model.
* **Fittings.** Users supply K sums; no fitting database is included.
* **ISO 7730 Table D.1:** three recalled rows could not be verified offline
  (above).
