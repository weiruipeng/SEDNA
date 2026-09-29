#!/usr/bin/env python3
"""
Filter Bamdam Results Tool

This script filters bamdam compute results and only keeps taxa that meet
specified criteria. It is designed to be integrated into the sedna-bamdam
workflow to filter taxa before plotting.

Usage:
    python filter_bamdam_results.py <tsv_file> <subs_file> <output_file> [--no-damage-filter] [--no-subs-filter] [--duplicity <float>] [--damage-min <float>] [--damage-max <float>] [--length <float>] [--ani <float>] [--uniq-kmer <float>]

Output:
    Saves a filtered TSV file with only taxa that pass all criteria.
    Also returns 0 (success) or 1 (failure) as exit code.
"""

import csv
import sys
import argparse
from io import StringIO


def parse_subs_file(subs_path):
    """Parse subs file to get 3' and 5' mutation rates at positions 1 and 5"""
    ga_3p_1 = None  # 3' end position 1 G->A
    ga_3p_5 = None  # 3' end position 5 G->A
    ct_5p_1 = None  # 5' end position 1 C->T (position -1 in bamdam output)
    ct_5p_5 = None  # 5' end position 5 C->T (position -5 in bamdam output)
    
    try:
        with open(subs_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Split by whitespace
                parts = line.split()
                if len(parts) < 2:
                    continue
                
                # Check for 3' end G->A at position 1 and 5
                if parts[0] == 'GA1':
                    try:
                        ga_3p_1 = float(parts[1])
                    except ValueError:
                        pass
                elif parts[0] == 'GA5':
                    try:
                        ga_3p_5 = float(parts[1])
                    except ValueError:
                        pass
                # Check for 5' end C->T at position 1 (position -1) and 5 (position -5)
                elif parts[0] == 'CT-1':
                    try:
                        ct_5p_1 = float(parts[1])
                    except ValueError:
                        pass
                elif parts[0] == 'CT-5':
                    try:
                        ct_5p_5 = float(parts[1])
                    except ValueError:
                        pass
    except FileNotFoundError:
        return None, None, None, None
    
    return ga_3p_1, ga_3p_5, ct_5p_1, ct_5p_5


def check_subs_filter(subs_path):
    """
    Check if the subs file passes the mutation rate filter.
    
    Returns:
        (passes, ga_3p_1, ga_3p_5, ct_5p_1, ct_5p_5)
    """
    ga_3p_1, ga_3p_5, ct_5p_1, ct_5p_5 = parse_subs_file(subs_path)
    
    # If any value is None, we can't apply the filter
    if ga_3p_1 is None or ga_3p_5 is None or ct_5p_1 is None or ct_5p_5 is None:
        return False, ga_3p_1, ga_3p_5, ct_5p_1, ct_5p_5
    
    # Check if 3' end position 1 > position 5 for G->A
    # and 5' end position 1 > position 5 for C->T
    passes = (ga_3p_1 > ga_3p_5) and (ct_5p_1 > ct_5p_5)
    
    return passes, ga_3p_1, ga_3p_5, ct_5p_1, ct_5p_5


def filter_tsv_file(tsv_path, subs_path, output_path, args):
    """
    Filter TSV file and write filtered results to output_path.
    
    Returns:
        (passed_taxa, total_taxa) - tuple of counts
    """
    total_taxa = 0
    passed_taxa = 0
    
    fieldnames_out = None
    rows_out = []
    
    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            total_taxa += 1
            passed = True
            ga_3p_1 = ga_3p_5 = ct_5p_1 = ct_5p_5 = None
            
            try:
                # Extract values from TSV
                duplicity = float(row['Duplicity'])
                damage_plus_1 = float(row['Damage+1'])
                damage_minus_1 = float(row['Damage-1'])
                mean_length = float(row['MeanLength'])
                ani = float(row['ANI'])
                uniq_kmers_per_read = float(row['UniqKmersPerRead'])
                
                # Check duplicity filter
                if duplicity >= args.duplicity:
                    passed = False
                
                # Check damage filters (if enabled)
                if not args.no_damage_filter:
                    if not (args.damage_min <= damage_plus_1 <= args.damage_max):
                        passed = False
                    if not (args.damage_min <= damage_minus_1 <= args.damage_max):
                        passed = False
                
                # Check length filter
                if mean_length >= args.length:
                    passed = False
                
                # Check ANI filter
                if ani <= args.ani:
                    passed = False
                
                # Check uniq kmer filter
                if uniq_kmers_per_read <= args.uniq_kmer:
                    passed = False
                
                # Check subs filter (if enabled)
                if passed and not args.no_subs_filter and subs_path:
                    subs_passes, ga_3p_1, ga_3p_5, ct_5p_1, ct_5p_5 = check_subs_filter(subs_path)
                    if not subs_passes:
                        passed = False
                
            except (ValueError, KeyError):
                passed = False
            
            # If passed, add to output
            if passed:
                passed_taxa += 1
                
                # Add subs data to row if available
                if ga_3p_1 is not None:
                    row['GA3p1'] = str(ga_3p_1)
                if ga_3p_5 is not None:
                    row['GA3p5'] = str(ga_3p_5)
                if ct_5p_1 is not None:
                    row['CT5p1'] = str(ct_5p_1)
                if ct_5p_5 is not None:
                    row['CT5p5'] = str(ct_5p_5)
                
                if fieldnames_out is None:
                    fieldnames_out = list(row.keys())
                
                rows_out.append(row)
    
    # Write output
    if fieldnames_out is not None:
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames_out, delimiter='\t')
            writer.writeheader()
            for row in rows_out:
                writer.writerow(row)
    
    return passed_taxa, total_taxa


def main():
    parser = argparse.ArgumentParser(
        description='Filter bamdam results based on specified criteria.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Default criteria:
  Duplicity < 2.0
  Damage+1 and Damage-1 between 0.1 and 0.6 (can be disabled with --no-damage-filter)
  MeanLength < 100.0
  ANI > 0.97
  UniqKmersPerRead > 0.4
  Subs filter: 3' position 1 > position 5 for G->A, 5' position 1 > position 5 for C->T (can be disabled with --no-subs-filter)
        """
    )
    
    parser.add_argument('tsv_file', help='Path to bamdam compute TSV file')
    parser.add_argument('subs_file', help='Path to bamdam compute subs file')
    parser.add_argument('output_file', help='Path to output filtered TSV file')
    
    # Filter thresholds
    parser.add_argument('--duplicity', type=float, default=2.0,
                        help='Maximum Duplicity (default: 2.0)')
    parser.add_argument('--damage-min', type=float, default=0.1,
                        help='Minimum Damage+1 and Damage-1 (default: 0.1)')
    parser.add_argument('--damage-max', type=float, default=0.6,
                        help='Maximum Damage+1 and Damage-1 (default: 0.6)')
    parser.add_argument('--length', type=float, default=100.0,
                        help='Maximum MeanLength (default: 100.0)')
    parser.add_argument('--ani', type=float, default=0.97,
                        help='Minimum ANI (default: 0.97)')
    parser.add_argument('--uniq-kmer', type=float, default=0.4,
                        help='Minimum UniqKmersPerRead (default: 0.4)')
    
    # Disable specific filters
    parser.add_argument('--no-damage-filter', action='store_true',
                        help='Disable damage-related filters (Damage+1, Damage-1)')
    parser.add_argument('--no-subs-filter', action='store_true',
                        help='Disable subs mutation rate filter')
    
    args = parser.parse_args()
    
    passed_taxa, total_taxa = filter_tsv_file(args.tsv_file, args.subs_file, 
                                            args.output_file, args)
    
    if passed_taxa == 0:
        print(f"ERROR: No taxa passed the filter out of {total_taxa} total taxa.")
        return 1
    else:
        print(f"Filtered {passed_taxa}/{total_taxa} taxa passed.")
        return 0


if __name__ == '__main__':
    sys.exit(main())
