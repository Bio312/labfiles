#!/usr/bin/env python3
"""
sum_pqr_charges.py — Sum net charges from PQR files robustly.

- Uses the second-to-last token on ATOM/HETATM lines as the charge
  (PQR ends with: x y z  charge  radius), so it tolerates variable columns.
- Handles .pqr and .pqr.gz
- Accepts files, directories, and globs; optional recursion.
- Outputs per-file results with columns: path, atoms, net_charge

New:
  - -o/--out path : write results to a file (CSV by default)
  - --sep csv|tsv : choose delimiter (default csv)
  - --no-header   : omit header row
"""

import argparse
import csv
import gzip
import os
import sys
from glob import glob

def is_pqr_line(s):
    return s.startswith("ATOM") or s.startswith("HETATM")

def iter_pqr_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    mode = "rt"
    try:
        with opener(path, mode, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                yield line.rstrip("\n")
    except OSError as e:
        sys.stderr.write("ERROR: cannot open {}: {}\n".format(path, e))

def sum_charge_for_file(path):
    """
    Returns (atom_count, net_charge_float).
    Charge = second-to-last token on ATOM/HETATM lines.
    """
    total = 0.0
    atoms = 0
    for line in iter_pqr_lines(path):
        if not is_pqr_line(line):
            continue
        parts = line.split()
        if len(parts) < 8:
            continue
        try:
            charge = float(parts[-2])
        except ValueError:
            continue
        atoms += 1
        total += charge
    return atoms, total

def collect_inputs(args):
    inputs = []
    for item in args.inputs:
        if os.path.isdir(item):
            if args.recursive:
                for root, _, files in os.walk(item):
                    for f in files:
                        if f.endswith(".pqr") or f.endswith(".pqr.gz"):
                            inputs.append(os.path.join(root, f))
            else:
                for f in os.listdir(item):
                    if f.endswith(".pqr") or f.endswith(".pqr.gz"):
                        inputs.append(os.path.join(item, f))
        else:
            matched = glob(item)
            if matched:
                inputs.extend([p for p in matched if os.path.isfile(p)])
            elif os.path.isfile(item):
                inputs.append(item)
            else:
                sys.stderr.write("WARNING: no such file or pattern: {}\n".format(item))
    # de-dup preserve order
    seen = set()
    unique = []
    for p in inputs:
        if p not in seen:
            unique.append(p)
            seen.add(p)
    return unique

def main():
    ap = argparse.ArgumentParser(description="Sum net charge from PQR files (robust to column variation).")
    ap.add_argument("inputs", nargs="+", help="PQR files, directories, or globs")
    ap.add_argument("-r", "--recursive", action="store_true", help="Recurse into directories")
    ap.add_argument("--warn-threshold", type=float, default=1000.0,
                    help="Absolute charge above which to warn (default: 1000)")
    ap.add_argument("--digits", type=int, default=3, help="Decimal places in output (default: 3)")
    ap.add_argument("-o", "--out", default="", help="Write results to this file (CSV by default)")
    ap.add_argument("--sep", choices=["csv", "tsv"], default="csv",
                    help="Output delimiter (csv=comma, tsv=tab). Default: csv")
    ap.add_argument("--no-header", action="store_true", help="Do not write header row")
    args = ap.parse_args()

    files = collect_inputs(args)
    if not files:
        sys.stderr.write("No input files found.\n")
        sys.exit(1)

    # formatter for net charge
    fmt = "{:." + str(args.digits) + "f}"

    # choose output stream + csv dialect
    delimiter = "," if args.sep == "csv" else "\t"
    outfh = open(args.out, "w", newline="", encoding="utf-8") if args.out else sys.stdout
    writer = csv.writer(outfh, delimiter=delimiter)

    # header
    if not args.no_header:
        writer.writerow(["path", "net_charge"])

    grand_atoms = 0
    grand_charge = 0.0

    for path in files:
        atoms, total = sum_charge_for_file(path)
        grand_atoms += atoms
        grand_charge += total

        if abs(total) > args.warn_threshold:
            sys.stderr.write(
                "WARNING: suspicious net charge ({}) in {}; this often means a parser assumed a fixed column index.\n"
                .format(total, os.path.basename(path))
            )

        writer.writerow([os.path.basename(path), fmt.format(total)])

    # Close file if we opened it (not stdout)
    if args.out:
        outfh.close()

    # Summary to stderr so it doesn't pollute the CSV/TSV
    sys.stderr.write("TOTAL\tatoms={}\tnet_charge={}\n".format(grand_atoms, fmt.format(grand_charge)))

if __name__ == "__main__":
    main()
