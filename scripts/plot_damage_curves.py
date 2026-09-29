#!/usr/bin/env python3
"""
Wrapper around `bamdam plotdamage` that generates standard postmortem damage
"smiley" plots for all taxa in a bamdam .subs file.

Usage:
  python plot_damage_curves.py --subs <subs_file> --tsv <tsv_file> --outdir <output_dir>
"""

import argparse
import csv
import os
import re
import subprocess
import sys


def load_taxa(tsv_path):
    """Load list of (TaxNodeID, TaxName) from bamdam TSV."""
    taxa = []
    with open(tsv_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            tid = row.get("TaxNodeID", "").strip()
            if not tid or tid in ("NA", "0"):
                continue
            name = row.get("TaxName", tid).strip().strip('"')
            taxa.append((tid, name))
    return taxa


def safe_filename(name):
    """Sanitize a taxon name for use in a filename."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:80] if s else "NA"


def main():
    ap = argparse.ArgumentParser(
        description="Generate damage plots for all taxa using bamdam plotdamage"
    )
    ap.add_argument("--subs", required=True, help="Path to bamdam .subs file")
    ap.add_argument("--tsv", required=True, help="Path to bamdam .tsv file (for taxon names)")
    ap.add_argument("--outdir", required=True, help="Output directory for plots")
    args = ap.parse_args()

    if not os.path.exists(args.subs):
        print(f"[plot_damage] .subs file not found: {args.subs}; skipping", file=sys.stderr)
        return

    # Verify bamdam is available
    try:
        subprocess.run(["bamdam", "plotdamage", "--help"], capture_output=True, timeout=10)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("[plot_damage] ERROR: bamdam plotdamage not found; skipping", file=sys.stderr)
        return

    taxa = load_taxa(args.tsv)
    if not taxa:
        print("[plot_damage] No taxa found in TSV; skipping", file=sys.stderr)
        return

    plot_dir = os.path.join(args.outdir, "damage_plots")
    os.makedirs(plot_dir, exist_ok=True)

    success = 0
    failed = 0
    for taxid, name in taxa:
        safe_name = safe_filename(name)
        out_png = os.path.join(plot_dir, f"{taxid}_{safe_name}_damage.png")

        # Skip if already exists
        if os.path.exists(out_png):
            success += 1
            continue

        cmd = [
            "bamdam", "plotdamage",
            "--in_subs", args.subs,
            "--tax", taxid,
            "--outplot", out_png,
        ]

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if result.returncode == 0 and os.path.exists(out_png):
                success += 1
            else:
                failed += 1
                print(f"[plot_damage] WARNING: bamdam plotdamage failed for {taxid} ({name}): {result.stderr.strip()}", file=sys.stderr)
        except subprocess.TimeoutExpired:
            failed += 1
            print(f"[plot_damage] WARNING: timeout for {taxid} ({name})", file=sys.stderr)
        except Exception as e:
            failed += 1
            print(f"[plot_damage] WARNING: error for {taxid} ({name}): {e}", file=sys.stderr)

    print(f"[plot_damage] Done: {success} plots generated, {failed} failed -> {plot_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
