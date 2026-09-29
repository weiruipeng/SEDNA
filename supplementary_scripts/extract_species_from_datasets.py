#!/usr/bin/env python3
import json
import sys
from pathlib import Path

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <datasets_unzip_dir> <output_species.txt>", file=sys.stderr)
    sys.exit(1)

root = Path(sys.argv[1])
outp = Path(sys.argv[2])

species = set()

TARGET_KEYS = {
    "organism_name",
    "organismName",
    "scientific_name",
    "scientificName",
    "sci_name",
    "sciName",
}

def walk(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k in TARGET_KEYS and isinstance(v, str):
                v = v.strip()
                if v:
                    species.add(v)
            walk(v)
    elif isinstance(obj, list):
        for x in obj:
            walk(x)

for fp in root.rglob("*"):
    if not fp.is_file():
        continue
    if fp.suffix.lower() not in {".json", ".jsonl"}:
        continue
    try:
        if fp.suffix.lower() == ".jsonl":
            with open(fp) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        obj = json.loads(line)
                        walk(obj)
                    except Exception:
                        pass
        else:
            with open(fp) as f:
                obj = json.load(f)
                walk(obj)
    except Exception:
        pass

with open(outp, "w") as fo:
    for sp in sorted(species):
        fo.write(sp + "\n")
