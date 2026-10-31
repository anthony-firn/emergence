import json, collections, sys
r = json.load(open(sys.argv[1]))
print({k: (len(v) if isinstance(v, list) else v) for k, v in r.items() if k not in ("source", "$schema", "date")})
c = collections.Counter((v["severity"], v["type"]) for v in r["violations"]); print(c.most_common(24))
c = collections.Counter(i["type"] for i in r.get("unconnected_items", [])); print("unconnected", c)
seen = collections.Counter()
for v in r["violations"] + r.get("unconnected_items", []):
    seen[v["type"]] += 1
    if seen[v["type"]] <= int(sys.argv[2]) : print(" ", v["type"], "|", v["description"][:96], "|", [(i["description"][:54], i["pos"]["x"], i["pos"]["y"]) for i in v["items"]][:2])
