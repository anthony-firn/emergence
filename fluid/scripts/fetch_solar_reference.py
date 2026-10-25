#!/usr/bin/env python3
"""Pinned public NSO/GONG maps; conservative 3x3 equal-area block averaging.

This reader supports only the archive's primary unscaled float32 image format.
It validates its WCS and units before converting to south-to-north SI maps.
No observation after the forecast origin is used in the initial map.
"""
from pathlib import Path
import gzip
import hashlib
import json
import math
import struct
import urllib.request

ROOT = Path(__file__).resolve().parents[1] / "crates/solar/data"
BASE = "https://gong2.nso.edu/oQR/zqs/"
FILES = [
    "201901/mrzqs190108/mrzqs190108t1114c2212_050.fits.gz",
    "201902/mrzqs190204/mrzqs190204t1104c2213_055.fits.gz",
]
HASHES = [
    "cace7a9d626f7e82c4cb368c358ead4f0dd6314f81acb5b0e6bfe8d2bd508e67",
    "fdef5ba15e3d8d3f90367acbfd34e409932b803100888ed4ede726d688320c69",
]


def prepare(path, expected_hash):
    url = BASE + path
    compressed = urllib.request.urlopen(url, timeout=30).read()
    if hashlib.sha256(compressed).hexdigest() != expected_hash:
        raise ValueError("GONG archive content changed: source hash does not match pinned observation")
    data = gzip.decompress(compressed)
    header = {}
    for offset in range(0, len(data), 80):
        card = data[offset:offset + 80].decode("ascii")
        if card.startswith("END "):
            image_start = ((offset + 80 + 2879) // 2880) * 2880
            break
        if card[8:10] == "= ":
            value = card[10:].split("/")[0].strip()
            header[card[:8].strip()] = value.strip("' ")
    else:
        raise ValueError("missing FITS END")
    assert int(header["BITPIX"]) == -32 and int(header["NAXIS"]) == 2
    assert int(header["NAXIS1"]) == 360 and int(header["NAXIS2"]) == 180
    assert header["CTYPE1"] == "CRLN-CEA" and header["CTYPE2"] == "CRLT-CEA"
    assert header["BUNIT"] == "Gauss"
    assert float(header.get("BSCALE", 1)) == 1 and float(header.get("BZERO", 0)) == 0
    assert float(header["CDELT1"]) == 1
    assert abs(float(header["CDELT2"]) - 2 / 180) < 1e-7
    assert float(header["CRPIX1"]) == 180.5 and float(header["CRPIX2"]) == 90.5
    assert float(header["CRVAL2"]) == 0
    initial_longitude = float(header["CRVAL1"]) - 180
    assert initial_longitude == int(initial_longitude)
    shift = int(initial_longitude) % 360
    values = struct.unpack_from(">" + "f" * (360 * 180), data, image_start)
    assert all(math.isfinite(v) for v in values)
    rows = []
    # Absolute Carrington longitudes 0.5,...,359.5 before averaging, rather
    # than treating the rolling archive map's first column as longitude zero.
    for latitude in range(0, 180, 3):
        row = []
        for longitude in range(0, 360, 3):
            total = sum(values[i * 360 + (j - shift) % 360]
                        for i in range(latitude, latitude + 3)
                        for j in range(longitude, longitude + 3))
            row.append(total / 9 * 1e-4)  # Gauss -> tesla
        rows.append(row)
    filename = "gong_" + header["DATE-OBS"].replace("-", "") + ".csv"
    text = "# 60x120 equal-area mu/longitude cell means, tesla, south to north\n"
    text += "\n".join(",".join(f"{v:.12e}" for v in row) for row in rows) + "\n"
    (ROOT / filename).write_text(text)
    return {
        "url": url,
        "compressed_sha256": hashlib.sha256(compressed).hexdigest(),
        "csv": filename,
        "csv_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "map_date_utc": header["DATE"],
        "original_shape": [180, 360], "prepared_shape": [60, 120],
        "original_first_longitude_edge_deg": shift,
        "units": "tesla",
    }


def main():
    ROOT.mkdir(parents=True, exist_ok=True)
    manifest = {
        "maps": [prepare(f, expected) for f, expected in zip(FILES, HASHES)],
        "preprocessing": "3x3 equal-area block means, Gauss to tesla; absolute Carrington longitude alignment",
        "scope": "Dated synoptic-map transport smoke comparison; maps are assembled over time, not simultaneous full-Sun measurements.",
        "instrument_notes": "GONG radial-field approximation and fitted polar fields; see NSO product documentation.",
        "source_documentation": "https://magmap.nso.edu/",
    }
    (ROOT / "gong_provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
