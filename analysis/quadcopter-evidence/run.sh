#!/usr/bin/env bash
# Re-runs every recorded software result for the quadcopter against the file as it stands.
# Needs Docker, a PX4 v1.17.0 checkout and a copy of AM32's Inc/targets.h. PX4 is built and
# flown with networking off and only its own folder mounted; see README.md.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../.." && pwd)
export EVIDENCE_WORK=${EVIDENCE_WORK:-/tmp/emergence-evidence}; W=$EVIDENCE_WORK; PX4=${PX4:?set PX4 to a PX4-Autopilot v1.17.0 checkout}
IMG=px4io/px4-dev:v1.16.0-rc1-258-g0369abd556; KICAD=kicad/kicad:9.0
mkdir -p "$W/kicad" "$W/px4-run"; cp "$HERE/fill.py" "$HERE/copper.py" "$W/kicad/"; cp "$HERE/fly.py" "$W/px4-run/"; chmod 777 "$W/kicad" "$W/px4-run"
if [ ! -x "$PX4/build/px4_sitl_nodds/bin/px4" ]; then
  printf 'CONFIG_MODULES_UXRCE_DDS_CLIENT=n\nCONFIG_MODULES_SIMULATION_GZ_BRIDGE=n\nCONFIG_MODULES_SIMULATION_GZ_MSGS=n\nCONFIG_MODULES_SIMULATION_GZ_PLUGINS=n\n' > "$PX4/boards/px4/sitl/nodds.px4board"
  docker run --rm --network none --memory 5g --cpus 4 --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PX4:/src" -w /src "$IMG" bash -c 'git config --global --add safe.directory "*"; cmake -S . -B build/px4_sitl_nodds -GNinja -DCONFIG=px4_sitl_nodds -DCMAKE_BUILD_TYPE=RelWithDebInfo && ninja -C build/px4_sitl_nodds -j4 px4 parameters_xml'
fi
cd "$ROOT/web"
node "$HERE/export-kicad.mjs" "$W/kicad/esc.kicad_pcb"
D="docker run --rm -v $W/kicad:/work -w /work --user $(id -u):$(id -g) -e HOME=/tmp"
$D --entrypoint python3 $KICAD /work/fill.py /work/esc.kicad_pcb /work/esc-filled.kicad_pcb
$D --entrypoint kicad-cli $KICAD pcb drc --format json --severity-all --all-track-errors -o /work/drc-filled.json /work/esc-filled.kicad_pcb
$D --entrypoint python3 $KICAD /work/copper.py /work/esc.kicad_pcb /work/copper.json
python3 "$HERE/drc-report.py" "$W/kicad/drc-filled.json" 2
node "$HERE/compare-copper.mjs"; node "$HERE/check-firmware-pins.mjs"; node "$HERE/px4-plant.mjs" "$W/px4-run/plant.json"
for tag in stock tuned; do docker run --rm --network none --memory 3g --cpus 3 --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PX4:/src" -v "$W/px4-run:/out" -w /src "$IMG" python3 /out/fly.py "$tag"; done
node "$HERE/px4-compare.mjs"; node "$HERE/record-evidence.mjs"
