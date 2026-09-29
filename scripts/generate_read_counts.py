#!/usr/bin/env python3
"""
Generate read-count summary table per taxonomic rank from bamdam TSV.

Usage:
  python generate_read_counts.py --tsv <bamdam.tsv> --rank <species|genus|family> --out <output.tsv>
"""

import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser(description="Generate read-count summary from bamdam TSV")
    ap.add_argument("--tsv", required=True, help="Path to bamdam .tsv file")
    ap.add_argument("--rank", required=True, help="Taxonomic rank label")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    rows = []
    with open(args.tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames_orig = reader.fieldnames or []
        for row in reader:
            taxid = row.get("TaxNodeID", "").strip()
            name = row.get("TaxName", "").strip().strip('"')
            total_reads = row.get("TotalReads", "0").strip()
            unagg_reads = row.get("UnaggregatedReads", "0").strip()
            rows.append({
                "TaxNodeID": taxid,
                "TaxName": name,
                "Rank": args.rank,
                "TotalReads": total_reads,
                "UnaggregatedReads": unagg_reads,
            })

    # Sort by TotalReads descending
    rows.sort(key=lambda r: int(r["TotalReads"]) if r["TotalReads"].isdigit() else 0, reverse=True)

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["TaxNodeID", "TaxName", "Rank", "TotalReads", "UnaggregatedReads"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    total_taxa = len(rows)
    total_reads = sum(int(r["TotalReads"]) for r in rows if r["TotalReads"].isdigit())
    print(f"[read_counts] {args.rank}: {total_taxa} taxa, {total_reads} total reads -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
