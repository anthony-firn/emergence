#!/usr/bin/env bash
# Project-local CAD converter: gmsh (with OpenCASCADE) in a virtualenv under
# external/tools/cad, used by geometry::cad to tessellate STEP/IGES/BREP files.
# Nothing is installed outside the repository.
set -euo pipefail
root="$(cd "$(dirname "$0")/../.." && pwd)"
venv="$root/external/tools/cad"
python3 -m venv "$venv"
"$venv/bin/pip" install --quiet --upgrade pip
"$venv/bin/pip" install --quiet "gmsh==4.13.1"
"$venv/bin/python" -c "import gmsh; gmsh.initialize(); print('gmsh', gmsh.__version__); gmsh.finalize()"
