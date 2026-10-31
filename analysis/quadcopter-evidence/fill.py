import sys, pcbnew
b = pcbnew.LoadBoard(sys.argv[1])
z = b.Zones(); print("zones", len(z), "kicad", pcbnew.Version())
pcbnew.ZONE_FILLER(b).Fill(z)
b.BuildConnectivity() if hasattr(b, "BuildConnectivity") else None
for x in z: print(x.GetNetname(), x.GetLayerName(), round(x.GetFilledArea() / 1e12, 1), "mm2", "filled" if x.IsFilled() else "EMPTY", "islands" , x.Outline().OutlineCount(), x.GetFilledPolysList(x.GetFirstLayer()).OutlineCount() if hasattr(x, "GetFilledPolysList") else "?")
pcbnew.SaveBoard(sys.argv[2], b)
