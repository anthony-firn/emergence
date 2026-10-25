#!/usr/bin/env python3
"""Build fluid/crates/compressible/data/elements_predicted.csv: estimated or
predicted properties of the elements whose melting and boiling points are
not both measured (mostly actinides and superheavy elements).

Source: the element infobox templates of the English Wikipedia
(Template:Infobox_<element>, raw wikitext), which quote the original
estimates and predictions with their references (relativistic DFT
calculations for Cn, Fl, Og; extrapolations for the late actinides and the
6d/7p elements). Every value is written with its qualifier as given there
('predicted', 'estimated', 'calculated', or none) and the reference names.

These values are model inputs of grade P (prediction): they are not data,
and nothing derived from them should be read as a measurement.
"""
import csv
import re
import sys
import urllib.request

UA = {"User-Agent": "Navier-Stokes-fluid-research/1.0 (element property builder)"}
OUT = "../crates/compressible/data/elements_predicted.csv"
BASE = "../crates/compressible/data/elements.csv"

FIELDS = {
    "t_melt": "melting point K",
    "t_boil": "boiling point K",
    "rho_solid": "density gpcm3nrt",
    "rho_liquid_mp": "density gpcm3mp",
    "h_fus": "heat fusion",
    "h_vap": "heat vaporization",
    "cp_molar": "heat capacity",
}


def fetch(title):
    url = "https://en.wikipedia.org/w/index.php?title=" + title + "&action=raw"
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
        return r.read().decode()


def field(text, key):
    """(value string, comment string) of infobox field `key`."""
    val = com = ""
    for line in text.split("\n"):
        m = re.match(r"\s*\|\s*" + re.escape(key) + r"\s*=(.*)$", line)
        if m:
            val = m.group(1).strip()
        # comment line: '<key> comment' or, for unit-suffixed keys such as
        # 'melting point K', '<key without the unit> comment'
        for ck in (key, key.rsplit(" ", 1)[0]):
            m = re.match(r"\s*\|\s*" + re.escape(ck) + r"\s+comment\s*=(.*)$", line)
            if m:
                com = m.group(1).strip()
    return val, com


def number(s):
    """First number in a wikitext value: '325 ± 15' -> 325, '8.91–8.95' ->
    8.93 (midpoint of a range), '{{nowrap|283 ± 11}}' -> 283."""
    s = re.sub(r"<ref[^>]*/>|<ref.*?</ref>", "", s)
    s = s.replace("{{nowrap|", "").replace("}}", "").replace("&nbsp;", " ").replace("−", "-")
    m = re.search(r"(-?\d+(?:\.\d+)?)\s*[–-]\s*(\d+(?:\.\d+)?)", s)
    if m and "±" not in s:
        return 0.5 * (float(m.group(1)) + float(m.group(2)))
    m = re.search(r"-?\d+(?:\.\d+)?", s)
    return float(m.group(0)) if m else None


def qualifier(com):
    c = com.lower()
    for q in ["predicted", "estimated", "extrapolated", "calculated"]:
        if q in c:
            return q
    return "quoted"


def refs(com):
    return ";".join(sorted(set(re.findall(r'<ref name\s*=\s*"?([^">/]+)"?', com))))


def main():
    base = list(csv.DictReader(open(BASE)))
    rows = []
    for el in base:
        if el["t_melt"] and el["t_boil"]:
            continue
        name = el["name"].lower()
        try:
            text = fetch("Template:Infobox_" + name)
        except Exception as e:  # noqa: BLE001
            print(f"{el['symbol']}: fetch failed ({e})", file=sys.stderr)
            continue
        row = {"z": el["z"], "symbol": el["symbol"]}
        notes = []
        for out_key, wiki_key in FIELDS.items():
            if el.get(out_key):
                continue  # measured value present: keep it
            val, com = field(text, wiki_key)
            x = number(val) if val else None
            if x is None:
                continue
            if out_key.startswith("rho"):
                x *= 1000.0  # g/cm³ -> kg/m³
            if out_key.startswith("h_"):
                x *= 1000.0  # kJ/mol -> J/mol (base table unit)
            row[out_key] = x
            notes.append(f"{out_key}:{qualifier(com)}[{refs(com)}]")
        if len(row) > 2:
            row["notes"] = " ".join(notes)
            rows.append(row)
            print(el["symbol"], row, file=sys.stderr)
    keys = ["z", "symbol"] + list(FIELDS) + ["notes"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    print(f"wrote {len(rows)} elements to {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
