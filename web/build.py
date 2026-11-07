"""Build the actual Rust models and package the static EvolveTech lab."""
import argparse
import hashlib
import json
import re
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def reusable_wasm(package):
    """Reuse an identified release only when its Rust sources are unchanged."""
    metadata = json.loads((package / "BUILD.json").read_text())
    commit = metadata.get("wasm_commit", metadata.get("commit", ""))
    wasm = package / "emergence_web.wasm"
    if not re.fullmatch(r"[a-f0-9]{40}", commit) or metadata.get("modified"):
        raise ValueError("A clean release with a full source commit is required.")
    if metadata.get("source", "").removesuffix(".git") != "https://github.com/anthony-firn/emergence":
        raise ValueError("The release must identify the Emergence source repository.")
    if hashlib.sha256(wasm.read_bytes()).hexdigest() != metadata.get("wasm_sha256"):
        raise ValueError("The release WASM checksum does not match BUILD.json.")
    subprocess.run(["git", "diff", "--exit-code", commit, "HEAD", "--", "fluid"], cwd=ROOT, check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", "fluid"], cwd=ROOT):
        raise ValueError("Rust sources have local changes; build the models instead.")
    return wasm, commit


def build(out, reuse_wasm_from=None):
    reused_commit = None
    if reuse_wasm_from:
        wasm, reused_commit = reusable_wasm(reuse_wasm_from.resolve())
    else:
        subprocess.run(["cargo", "build", "--locked", "--profile", "web", "-p", "emergence-web", "--target", "wasm32-unknown-unknown"], cwd=ROOT / "fluid", check=True)
        wasm = ROOT / "fluid/target/wasm32-unknown-unknown/web/emergence_web.wasm"
    out.mkdir(parents=True, exist_ok=True)
    for name in ("index.html", "lab.css", "lab.js", "worker.js"):
        shutil.copy2(ROOT / "web" / name, out / name)
    # The Node agent and unused vendored builds are not part of the site.
    shutil.copytree(ROOT / "web" / "cad", out / "cad", dirs_exist_ok=True, ignore=shutil.ignore_patterns("agent", "vendor", "__pycache__"))
    shutil.copy2(wasm, out / wasm.name)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT))
    provenance = {"source": "https://github.com/anthony-firn/emergence", "commit": commit,
                  "modified": dirty, "wasm_sha256": hashlib.sha256(wasm.read_bytes()).hexdigest()}
    if reused_commit:
        provenance["wasm_commit"] = reused_commit
    (out / "BUILD.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"Packaged emergence lab into {out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "web/dist")
    parser.add_argument("--reuse-wasm-from", type=Path, help="Verified prior package (BUILD.json + WASM); rejects changed Rust sources")
    args = parser.parse_args()
    build(args.out.resolve(), args.reuse_wasm_from)
