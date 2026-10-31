# Quadcopter evidence

The scripts behind the *Real software* results in the quadcopter's test plan. `run.sh` repeats
them against `web/cad/projects/inspection-drone.json` and records the results through the
project's own commands, so each one carries a digest of the documents it covers. Change any
of those documents and the result reads "stale" until this is run again.

| Check | Tool | What it does |
|---|---|---|
| Parameters | PX4 v1.17.0 `parameters.json` | Every parameter the project uses exists, with the same default and within range |
| Flight | PX4 v1.17.0 in SIH | Flies this aircraft's mass, inertia, thrust slope and motor lag, with a roll step, 24 averaged rate steps and a frequency sweep, and compares the log with the project's own simulator on the same plant |
| Board | KiCad 9.0.9 | Fills zones with KiCad's filler and runs its design-rule check on the exported board |
| Copper | KiCad 9.0.9 | Compares the app's Gerbers with the copper KiCad computes, per layer, on a 50 µm grid |
| Pins | AM32 source, read as text | Compares every pin AM32 drives with the net on the board, using pin numbers from KiCad's symbol |

**Not yet done:** `am32-target.h` is this board's firmware target. It has been written and its
pins checked as text, but it has not been compiled. Until it is, `am32-build` reads "not run".

**SIH is not the aircraft.** Its thrust is linear in command, so `SIH_T_MAX` is set to the
real aircraft's slope at hover to keep the same control authority, and the project's simulator
is run on that same plant for the comparison. It injects heavy gyro noise, which is why rate
steps are averaged. Its control loop runs at 250 Hz, where a Pixhawk 6X runs at 400 Hz.

**Isolation.** PX4 is compiled and run in its own published image with `--network none` and only
its own folder mounted. Its DDS and Gazebo modules are switched off because they fetch code at
build time.

`esc-board.mjs` is the layout as commands, kept for how the board was made. The project file is
the source of truth: `domains/rebuild.js` writes any project back out as commands from blank.
