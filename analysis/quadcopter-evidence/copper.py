import sys, json, pcbnew
b = pcbnew.LoadBoard(sys.argv[1]); pcbnew.ZONE_FILLER(b).Fill(b.Zones())
print("zone layers", sorted({pcbnew.LayerName(l) for z in b.Zones() for l in z.GetLayerSet().Seq()}), [ (z.GetNetname(), pcbnew.LayerName(z.GetLayerSet().Seq()[0])) for z in b.Zones()][:2])
out = {}
for name, layer in (("top", pcbnew.F_Cu), ("bottom", pcbnew.B_Cu)):
    s = pcbnew.SHAPE_POLY_SET(); b.ConvertBrdLayerToPolygonalContours(layer, s); s.Simplify() if hasattr(s, "Simplify") else None
    polys = []
    for i in range(s.OutlineCount()):
        o = s.Outline(i); rings = [[(o.CPoint(k).x / 1e6, o.CPoint(k).y / 1e6) for k in range(o.PointCount())]]
        for h in range(s.HoleCount(i)):
            r = s.Hole(i, h); rings.append([(r.CPoint(k).x / 1e6, r.CPoint(k).y / 1e6) for k in range(r.PointCount())])
        polys.append(rings)
    out[name] = polys; print(name, len(polys), "polygons", sum(len(p) - 1 for p in polys), "holes")
json.dump(out, open(sys.argv[2], "w"))
