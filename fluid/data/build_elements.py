#!/usr/bin/env python3
"""Build fluid/crates/compressible/data/elements.csv: per-element properties
for the phase models, every value tagged with its source.

Sources (fetched at build time):
  PubChem periodic table        https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV
  mendeleev (MIT) database      pip install mendeleev; values from the CRC
                                Handbook (Haynes 2014/2016) with citation keys
  Wikipedia element data pages  (raw wikitext; values from CRC, WebElements,
                                Kaye & Laby, as quoted there):
     Elastic_properties_of_the_elements_(data_page)
     Speeds_of_sound_of_the_elements
     Densities_of_the_elements_(data_page)

Usage: python3 build_elements.py OUT.csv   (needs network and mendeleev)
Units in the output are SI: kg/m³, K, J/mol, J/(mol K), Pa, m/s.
"""
import csv
import io
import re
import sys
import urllib.request

UA = {"User-Agent": "fluid-sim-research/0.1 (element property compilation)"}


def fetch(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8")


def wiki(page):
    return fetch(f"https://en.wikipedia.org/w/index.php?title={page}&action=raw")


NUM = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?")


def first_number(s):
    s = s.replace("&nbsp;", " ").replace("−", "-").replace(",", "")
    s = re.sub(r"<ref[^>]*>.*?</ref>|<ref[^/]*/>", "", s)
    m = NUM.search(s)
    return float(m.group()) if m else None


def parse_elastic(text):
    """{section: {Z: value}} for the four elastic tables (GPa or ratio)."""
    out = {}
    sections = re.split(r"^== (.+?) ==\s*$", text, flags=re.M)
    for title, body in zip(sections[1::2], sections[2::2]):
        vals = {}
        lines = body.splitlines()
        for i, line in enumerate(lines):
            m = re.match(r"^\|\s*(\d+)\s*\|\|\s*([A-Z][a-z]?)\s*\|\|", line)
            if m and i + 1 < len(lines):
                use = lines[i + 1].lstrip("|").split("||")[0]
                v = first_number(use)
                if v is not None:
                    vals[int(m.group(1))] = v
        out[title.strip()] = vals
    return out


def blocks(text, header_re):
    """Yield (Z, label, [rows]) for '! colspan... <!-- Z Sym --> ...' blocks."""
    cur = None
    for line in text.splitlines():
        m = re.match(header_re, line)
        if m:
            if cur:
                yield cur
            cur = (int(m.group(1)), line, [])
            continue
        if cur and line.startswith("|") and not line.startswith("|-") and not line.startswith("|}"):
            cur[2].append(line)
    if cur:
        yield cur


def parse_sound(text):
    solid, liquid = {}, {}
    solid_part, _, fluid_part = text.partition("==Speed of sound, fluid phases==")
    hdr = r"^!\s*colspan=\"\d\"\s*\|\s*<!--\s*(\d+)\s*[A-Za-z]+\s*-->"
    for z, label, rows in blocks(solid_part, hdr):
        for r in rows:
            cells = [c.strip() for c in r.lstrip("|").split("||")]
            if cells[0] == "use":
                lon = first_number(cells[1]) if len(cells) > 1 else None
                tr = first_number(cells[2]) if len(cells) > 2 else None
                ext = first_number(cells[3]) if len(cells) > 3 else None
                solid[z] = (lon, tr, ext)
    for z, label, rows in blocks(fluid_part, hdr):
        if "(liquid)" not in label:
            continue
        for r in rows:
            cells = [c.strip() for c in r.lstrip("|").split("||")]
            if cells[0] == "use" and len(cells) > 1:
                v = first_number(cells[1])
                note = cells[2] if len(cells) > 2 else ""
                tc = first_number(note) if "°C" in note else None
                if v is not None:
                    liquid[z] = (v, tc + 273.15 if tc is not None else None)
    return solid, liquid


def density_value(s):
    """The number written directly before its unit (notes such as
    '(at −35 °C)' must not be read as the value)."""
    s = s.replace("&nbsp;", " ")
    m = re.search(r"(\d+(?:\.\d+)?)\s*(g/cm|g/L|g/l|kg/m)", s)
    if not m:
        return None
    v = float(m.group(1))
    return v * 1000.0 if m.group(2) == "g/cm" else v


def parse_use_blocks(text):
    """{Z: value} from '! colspan="2" | <!-- Z Sym -->' blocks' 'use' rows."""
    hdr = r"^!\s*colspan=\"2\"[^|]*\|\s*<!--\s*(\d+)\s*[A-Za-z]+\s*-->"
    out = {}
    for z, label, rows in blocks(text, hdr):
        for r in rows:
            cells = [c.strip() for c in r.lstrip("|").split("||")]
            if cells[0] == "use" and len(cells) > 1:
                v = first_number(re.sub(r"\([^)]*\)", "", cells[1]))
                if v is not None and z not in out:
                    out[z] = v
    return out


def parse_density(text):
    sol, _, rest = text.partition("== Density, liquid phase ==")
    liq, _, _ = rest.partition("== Density, gas phase ==")
    hdr = r"^!\s*colspan=\"2\"[^|]*\|\s*<!--\s*(\d+)\s*[A-Za-z]+\s*-->"
    out = []
    for part in (sol, liq):
        d = {}
        for z, label, rows in blocks(part, hdr):
            for r in rows:
                cells = [c.strip() for c in r.lstrip("|").split("||")]
                if cells[0] == "use" and len(cells) > 1:
                    v = density_value(cells[1])
                    if v is not None:
                        d[z] = v
        out.append(d)
    return out


def parse_vapor_pressure(text):
    """{Z: [(p_Pa, T_K)]} from the 'use' rows of the vapour-pressure table."""
    head = text.split("|-", 2)[0] if False else text
    pressures = []
    for line in text.splitlines():
        m = re.match(r"^!\s*rowspan=\"2\"[^|]*\|\s*([\d.]+)\s*(k|M)?\s*$", line)
        if m:
            v = float(m.group(1)) * {None: 1.0, "k": 1e3, "M": 1e6}[m.group(2)]
            pressures.append(v)
        if line.startswith("! colspan"):
            break
    out = {}
    hdr = r"^!\s*colspan=\"\d+\"\s*\|\s*<!--\s*(\d+)\s*[A-Za-z]+\s*-->"
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = re.match(hdr, lines[i])
        if m:
            z = int(m.group(1))
            j = i + 1
            while j < len(lines) and not lines[j].startswith("| use"):
                if re.match(hdr, lines[j]):
                    break
                j += 1
            if j < len(lines) and lines[j].startswith("| use"):
                vals = []
                for k in range(len(pressures)):
                    if j + 1 + k >= len(lines):
                        break
                    cell = re.sub(r"<!--.*?-->", "", lines[j + 1 + k]).lstrip("|")
                    v = first_number(cell)
                    if v is not None:
                        vals.append((pressures[k], v))
                if vals and z not in out:
                    out[z] = vals
            i = j
        i += 1
    return out, pressures


INFOBOX_NAME = {"Aluminum": "aluminium", "Cesium": "caesium"}


def infobox_densities(name):
    """Liquid density at m.p. (or else b.p.) and solid density near room
    temperature from Template:Infobox <element>, kg/m³."""
    import time
    time.sleep(0.3)
    try:
        t = wiki("Template:Infobox_" + INFOBOX_NAME.get(name, name.lower()))
    except Exception:
        return None, None, None
    def field(key):
        m = re.search(r"^\|" + re.escape(key) + r"\s*=\s*([^\n|]*)", t, flags=re.M)
        return first_number(m.group(1)) if m else None
    mp, bp, rt = field("density gpcm3mp"), field("density gpcm3bp"), field("density gpcm3nrt")
    liq = mp if mp is not None else bp
    tag = "infobox(mp)" if mp is not None else ("infobox(bp)" if bp is not None else None)
    return (liq * 1000.0 if liq is not None else None), (rt * 1000.0 if rt is not None else None), tag


def main(out_path):
    pub = list(csv.DictReader(io.StringIO(fetch("https://pubchem.ncbi.nlm.nih.gov/rest/pug/periodictable/CSV"))))
    from mendeleev.fetch import fetch_table

    el = fetch_table("elements").set_index("atomic_number")
    pt = fetch_table("phasetransitions")
    elastic = parse_elastic(wiki("Elastic_properties_of_the_elements_(data_page)"))
    s_solid, s_liquid = parse_sound(wiki("Speeds_of_sound_of_the_elements"))
    rho_s_wiki, rho_l_wiki = parse_density(wiki("Densities_of_the_elements_(data_page)"))
    hfus_wiki = parse_use_blocks(wiki("Heats_of_fusion_of_the_elements"))
    hvap_wiki = parse_use_blocks(wiki("Heats_of_vaporization_of_the_elements_(data_page)"))
    alin_wiki = parse_use_blocks(wiki("Thermal_expansivities_of_the_elements"))

    def pick_transition(z):
        rows = pt[(pt.atomic_number == z) & (pt.is_transition != 1.0)]
        if rows.empty:
            return None
        # prefer the row without an allotrope label, else the first listed
        plain = rows[rows.allotrope.isna()]
        return (plain if not plain.empty else rows).iloc[0]

    cols = [
        "z", "symbol", "name", "mass_u", "standard_state", "block", "radioactive",
        "rho_solid", "rho_liquid_mp", "t_melt", "t_boil", "sublimes", "allotrope",
        "h_fus", "h_vap", "cp_molar", "t_crit", "p_crit", "t_triple", "p_triple",
        "bulk_modulus", "shear_modulus", "c_long_solid", "c_trans_solid",
        "c_liquid", "c_liquid_t", "alpha_lin_solid", "data_flags", "sources",
    ]
    rows = []
    for p in pub:
        z = int(p["AtomicNumber"])
        e = el.loc[z]
        src = []

        def num(x):
            try:
                v = float(x)
                return v if v == v else None
            except (TypeError, ValueError):
                return None

        t = pick_transition(z)
        t_melt = num(t.melting_point) if t is not None else None
        t_boil = num(t.boiling_point) if t is not None else None
        if t_melt is None:
            t_melt = num(p["MeltingPoint"])
            if t_melt is not None:
                src.append("t_melt:pubchem")
        else:
            src.append("t_melt:haynes2016crc")
        if t_boil is None:
            t_boil = num(p["BoilingPoint"])
            if t_boil is not None:
                src.append("t_boil:pubchem")
        else:
            src.append("t_boil:haynes2016crc")
        sub = bool(t is not None and t.is_sublimation_point == 1.0)
        rho_s = rho_s_wiki.get(z)
        if rho_s is not None:
            src.append("rho_solid:wiki-densities(use)")
        elif num(e.density) is not None and "Gas" not in p["StandardState"]:
            # for elements that are gases at room temperature mendeleev's
            # density is the gas density, not the solid's: leave it empty
            rho_s = num(e.density) * 1000.0
            src.append("rho_solid:mendeleev(haynes2014crc)")
        rho_l = rho_l_wiki.get(z)
        if rho_l is not None:
            src.append("rho_liquid_mp:wiki-densities(use)")
        elif t_melt is not None:
            li, so, tag = infobox_densities(p["Name"])
            if li is not None:
                rho_l = li
                src.append(f"rho_liquid:{tag}")
            if rho_s is None and so is not None and "Gas" not in p["StandardState"]:
                rho_s = so
                src.append("rho_solid:infobox(nrt)")
        flags = []

        def latent(wv, mv, key):
            # CRC values quoted on Wikipedia first; mendeleev as fallback;
            # flag disagreements (mendeleev had K's two heats swapped)
            if wv is not None:
                src.append(f"{key}:wiki-heats(use)")
                if mv is not None and abs(mv / wv - 1.0) > 0.2:
                    flags.append(f"{key}: mendeleev {mv:g} vs CRC {wv:g} kJ/mol")
                return wv
            if mv is not None:
                src.append(f"{key}:mendeleev")
            return mv

        hf = latent(hfus_wiki.get(z), num(e.fusion_heat), "h_fus")
        hv = latent(hvap_wiki.get(z), num(e.evaporation_heat), "h_vap")
        cp = num(e.molar_heat_capacity)
        alin = alin_wiki.get(z)
        if alin is not None:
            src.append("alpha_lin:wiki-expansivities(use)")
        if cp is not None:
            src.append("cp_molar:mendeleev(haynes2014crc)")
        K = elastic.get("Bulk modulus", {}).get(z)
        G = elastic.get("Shear modulus", {}).get(z)
        if K is not None:
            src.append("bulk_modulus:wiki-elastic(use)")
        if G is not None:
            src.append("shear_modulus:wiki-elastic(use)")
        cl, ct, _ = s_solid.get(z, (None, None, None))
        if cl or ct:
            src.append("c_solid:wiki-sound(use)")
        liq = s_liquid.get(z)
        if liq:
            src.append("c_liquid:wiki-sound(use)")

        def f(v, scale=1.0):
            return "" if v is None else f"{v * scale:.6g}"

        rows.append({
            "z": z, "symbol": p["Symbol"], "name": p["Name"], "mass_u": p["AtomicMass"],
            "standard_state": p["StandardState"], "block": p["GroupBlock"],
            "radioactive": int(bool(e.is_radioactive)),
            "rho_solid": f(rho_s), "rho_liquid_mp": f(rho_l),
            "t_melt": f(t_melt), "t_boil": f(t_boil), "sublimes": int(sub),
            "allotrope": "" if t is None or not isinstance(t.allotrope, str) else t.allotrope,
            "h_fus": f(hf, 1000.0), "h_vap": f(hv, 1000.0), "cp_molar": f(cp),
            "t_crit": f(num(t.critical_temperature) if t is not None else None),
            "p_crit": f(num(t.critical_pressure) if t is not None else None, 1e6),
            "t_triple": f(num(t.triple_point_temperature) if t is not None else None),
            "p_triple": f(num(t.triple_point_pressure) if t is not None else None, 1e3),
            "bulk_modulus": f(K, 1e9), "shear_modulus": f(G, 1e9),
            "c_long_solid": f(cl), "c_trans_solid": f(ct),
            "c_liquid": f(liq[0] if liq else None), "c_liquid_t": f(liq[1] if liq else None),
            "alpha_lin_solid": f(alin, 1e-6), "data_flags": " | ".join(flags),
            "sources": ";".join(src),
        })
    vp, pressures = parse_vapor_pressure(wiki("Vapor_pressures_of_the_elements"))
    vp_path = out_path.replace("elements.csv", "vapor_pressure.csv")
    with open(vp_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["z", "p_pa", "t_k"])
        for z in sorted(vp):
            for p_, t_ in vp[z]:
                w.writerow([z, f"{p_:g}", f"{t_:g}"])
    print(f"vapour pressure: {len(vp)} elements, pressure levels {pressures}")
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    have = lambda k: sum(1 for r in rows if r[k] != "")
    print(f"wrote {len(rows)} elements to {out_path}")
    for r in rows:
        if r["data_flags"]:
            print(f"  flag {r['symbol']}: {r['data_flags']}")
    for k in ["alpha_lin_solid", "rho_solid", "rho_liquid_mp", "t_melt", "t_boil", "h_fus", "h_vap", "cp_molar", "t_crit", "t_triple", "bulk_modulus", "shear_modulus", "c_liquid"]:
        print(f"  {k:15} {have(k):3} / {len(rows)}")


if __name__ == "__main__":
    main(sys.argv[1])
