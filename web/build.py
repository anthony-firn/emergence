"""Build the actual Rust models and package the static EvolveTech lab."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def build(out):
    subprocess.run(["cargo", "build", "--locked", "--profile", "web", "-p", "emergence-web", "--target", "wasm32-unknown-unknown"], cwd=ROOT / "fluid", check=True)
    out.mkdir(parents=True, exist_ok=True)
    for name in ("index.html", "lab.css", "lab.js", "worker.js"):
        shutil.copy2(ROOT / "web" / name, out / name)
    shutil.copytree(ROOT / "web" / "cad", out / "cad", dirs_exist_ok=True)
    wasm = ROOT / "fluid/target/wasm32-unknown-unknown/web/emergence_web.wasm"
    shutil.copy2(wasm, out / wasm.name)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT))
    provenance = {"source": "https://github.com/anthony-firn/emergence", "commit": commit,
                  "modified": dirty, "wasm_sha256": hashlib.sha256(wasm.read_bytes()).hexdigest()}
    (out / "BUILD.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Packaged emergence lab into {out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "web/dist")
    build(parser.parse_args().out.resolve())
