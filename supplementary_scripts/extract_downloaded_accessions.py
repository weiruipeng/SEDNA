#!/usr/bin/env python3
import json
import os
import re
import sys
from pathlib import Path

# 兼容所有Python版本，去掉高版本类型注解
ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")

def extract_from_text_file(fp, results):
    try:
        text = fp.read_text(errors="ignore")
    except Exception:
        return
    for m in ACC_RE.finditer(text):
        results.add(m.group(1))

def extract_from_json_file(fp, results):
    try:
        data = json.loads(fp.read_text(errors="ignore"))
    except Exception:
        return

    def walk(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if isinstance(v, (dict, list)):
                    walk(v)
                elif isinstance(v, str):
                    for m in ACC_RE.finditer(v):
                        results.add(m.group(1))
        elif isinstance(obj, list):
            for item in obj:
                walk(item)
        elif isinstance(obj, str):
            for m in ACC_RE.finditer(obj):
                results.add(m.group(1))

    walk(data)

def main():
    if len(sys.argv) != 3:
        print("Usage: {} <target_genomes_dir> <output_acc_list>".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)

    target_dir = Path(sys.argv[1])
    output_file = Path(sys.argv[2])

    if not target_dir.exists():
        print("[ERROR] directory not found: {}".format(target_dir), file=sys.stderr)
        sys.exit(1)

    results = set()

    for root, _, files in os.walk(target_dir):
        for fn in files:
            fp = Path(root) / fn
            if ACC_RE.search(fn):
                for m in ACC_RE.finditer(fn):
                    results.add(m.group(1))

            lower = fn.lower()
            if lower.endswith(".json"):
                extract_from_json_file(fp, results)
            elif lower.endswith((".txt", ".tsv", ".csv", ".fna", ".fa", ".fasta", ".gbff", ".gff", ".md")):
                extract_from_text_file(fp, results)

    with output_file.open("w") as out:
        for acc in sorted(results):
            out.write(acc + "\n")

if __name__ == "__main__":
    main()
