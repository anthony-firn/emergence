"""Export a labelled scientific preview and movies from the Rust renderer.

python3 analysis/export_solar_render.py /tmp/solar-render-fixed analysis/solar_render /tmp/solar-render-orbit
Requires matplotlib and ffmpeg; they are not simulation dependencies.
"""
from pathlib import Path
import json
import shutil
import subprocess
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def movie(source, target, orbit):
    font = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
    metadata = json.loads((source / "metadata.json").read_text())
    step_days = float(metadata["days"]) / max(int(metadata["frames"]) - 1, 1)
    title = "Inspection camera orbit; not solar rotation" if orbit else "Fixed Carrington view; simulated surface transport"
    filters = (
        "[0:v][1:v]hstack=inputs=2:shortest=1,pad=iw:1100:0:150:color=0x080c12,"
        f"drawtext=fontfile={font}:text='SOLAR CONTINUUM AND MAGNETIC FIELD':x=60:y=24:fontsize=36:fontcolor=white,"
        f"drawtext=fontfile={font}:text='{title}':x=60:y=72:fontsize=24:fontcolor=0xaabccc,"
        f"drawtext=fontfile={font}:text='Grey LTE continuum - close view':x=60:y=116:fontsize=20:fontcolor=white,"
        f"drawtext=fontfile={font}:text='PFSS geometry - wider view - no coronal emission':x=960:y=116:fontsize=20:fontcolor=white,"
        f"drawtext=fontfile={font}:text='Day %{{expr\\:floor(t*4)*{step_days:.12g}}} - Gold closed / green open - no new magnetic emergence':x=60:y=1060:fontsize=22:fontcolor=white,"
        # Retain the last source snapshot for its full quarter-second. Without
        # padding, frame-rate conversion ends at the last timestamp and drops
        # the final simulated day's display interval.
        "tpad=stop_mode=clone:stop_duration=0.25,fps=24[out]"
    )
    subprocess.run([
        "ffmpeg", "-hide_banner", "-loglevel", "error", "-y",
        "-filter_complex_threads", "4", "-framerate", "4", "-loop", "1",
        "-i", str(source / "continuum.png"), "-framerate", "4",
        "-i", str(source / "frame%04d.png"), "-filter_complex", filters,
        "-map", "[out]", "-c:v", "libx264", "-threads", "8", "-preset", "medium",
        "-crf", "19", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(target)
    ], check=True)


def main():
    source = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((source / "metadata.json").read_text())
    fig = plt.figure(figsize=(14, 8), facecolor="#080c12")
    images = ["continuum.png", "frame0000.png"]
    titles = ["Visible continuum · grey LTE atmosphere", "Radial field and potential coronal topology"]
    for x, name, title in zip([0.025, 0.525], images, titles):
        ax = fig.add_axes([x, 0.17, 0.45, 0.69])
        ax.imshow(plt.imread(source / name))
        ax.axis("off")
        ax.set_title(title, color="#dce7f2", fontsize=12, pad=8)
    fig.text(0.04, 0.95, "SOLAR CONTINUUM AND MAGNETIC FIELD", color="white", fontsize=21)
    fig.text(0.04, 0.905, "GONG initial boundary · 8 January 2019 · Day 0 · different view scales", color="#aabccc", fontsize=12)
    fig.text(0.04, 0.125, f"Continuum from spectral transfer · effective temperature {metadata['effective_temperature_K']:.1f} K", color="#aabccc", fontsize=11)
    fig.text(0.53, 0.125, "Surface field: blue negative / red positive · asinh scale ±20 G", color="#aabccc", fontsize=11)
    fig.text(0.53, 0.08, "Gold: closed curves      Green: open curves", color="#f3c477", fontsize=11)
    fig.text(0.04, 0.03, "Magnetic curves show geometry, not coronal emission. No new emergence, spots, heating, reconnection, or eruption forecast.", color="#aabccc", fontsize=10)
    fig.savefig(output / "preview.png", dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)
    for name in ["metadata.json", "frames.csv", "continuum.png", "frame0000.png", "frame0027.png"]:
        if (source / name).exists():
            shutil.copyfile(source / name, output / name)
    movie(source, output / "transport.mp4", False)
    if len(sys.argv) > 3:
        orbit = Path(sys.argv[3]).resolve()
        movie(orbit, output / "inspection.mp4", True)
        shutil.copyfile(orbit / "metadata.json", output / "inspection_metadata.json")
    viewer = """<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Solar render</title>
<style>body{max-width:1400px;margin:30px auto;padding:20px;background:#080c12;color:#dce7f2;font:17px system-ui}video{width:100%}button{padding:10px 18px;margin:10px 10px 10px 0;background:#19313e;color:white;border:1px solid #476777;border-radius:6px}p{line-height:1.6}a{color:#8cd8c5}</style>
<h1>Solar continuum and magnetic field</h1>
<p>GONG boundary · 8 January 2019 · 27 days of prescribed surface transport without new emergence.</p>
<button onclick="choose('transport.mp4','Fixed Carrington camera')">Surface transport</button>
<button id="orbit" onclick="choose('inspection.mp4','Inspection camera orbit; not solar rotation')">Inspection camera</button>
<p id="caption">Fixed Carrington camera</p>
<video id="movie" controls preload="metadata" poster="preview.png" src="transport.mp4"></video>
<p>The continuum is a grey LTE transfer calculation. Gold/green magnetic curves show closed/open geometry, not coronal emission. The views have different scales. The potential corona has no currents, heating or reconnection; these movies do not predict eruptions.</p>
<p><a href="metadata.json">Provenance and model settings</a> · <a href="frames.csv">Per-frame diagnostics</a></p>
<script>const movie=document.getElementById('movie');function choose(path,caption){movie.src=path;document.getElementById('caption').textContent=caption;movie.load()}__ORBIT__</script></html>"""
    viewer = viewer.replace("__ORBIT__", "" if len(sys.argv) > 3 else "document.getElementById('orbit').hidden=true;")
    (output / "index.html").write_text(viewer)
    print(output)


if __name__ == "__main__":
    main()
