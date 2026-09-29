#!/usr/bin/env python3
import gzip
import sys
from pathlib import Path

if len(sys.argv) < 4:
    print(
        f"Usage: {sys.argv[0]} <db_accessions.txt> <output.txt> <acc2taxid.gz ...>",
        file=sys.stderr
    )
    sys.exit(1)

acc_file = Path(sys.argv[1])
out_file = Path(sys.argv[2])
map_files = [Path(x) for x in sys.argv[3:]]

wanted = {}
with open(acc_file) as f:
    for line in f:
        acc = line.strip()
        if not acc:
            continue
        wanted[acc] = None
        short = acc.split(".")[0]
        if short not in wanted:
            wanted[short] = None

resolved = {}

for mf in map_files:
    with gzip.open(mf, "rt") as f:
        header = next(f, None)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc = parts[0]          # accession
            acc_ver = parts[1]      # accession.version
            taxid = parts[2]

            if acc_ver in wanted and acc_ver not in resolved:
                resolved[acc_ver] = taxid
            if acc in wanted and acc not in resolved:
                resolved[acc] = taxid

# 输出时按原数据库 accessions 列表顺序来
with open(acc_file) as fi, open(out_file, "w") as fo:
    for line in fi:
        orig = line.strip()
        if not orig:
            continue
        taxid = resolved.get(orig)
        if taxid is None:
            taxid = resolved.get(orig.split(".")[0])
        if taxid is not None:
            fo.write(f"{orig}\t{orig}\t{taxid}\t0\n")
