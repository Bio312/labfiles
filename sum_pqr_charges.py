#!/usr/bin/env python3
import os
import sys
import glob

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python3 sum_pqr_charges.py <input_dir_with_pqr> <output_tsv>\n")
        sys.exit(1)

    in_dir  = os.path.expanduser(os.path.expandvars(sys.argv[1]))
    out_tsv = os.path.expanduser(os.path.expandvars(sys.argv[2]))

    if not os.path.isdir(in_dir):
        sys.stderr.write(f"ERROR: Not a directory: {in_dir}\n")
        sys.exit(1)

    pqr_files = sorted(glob.glob(os.path.join(in_dir, "*.pqr")))
    if not pqr_files:
        sys.stderr.write(f"ERROR: No .pqr files found in {in_dir}\n")
        sys.exit(1)

    with open(out_tsv, "w") as outfh:
        for f in pqr_files:
            total = 0.0
            with open(f) as fh:
                for line in fh:
                    if line.startswith(("ATOM", "HETATM")):
                        parts = line.split()
                        # Column 10 (index 9) = charge
                        if len(parts) >= 10:
                            try:
                                total += float(parts[9])
                            except ValueError:
                                continue
            outfh.write(f"{os.path.basename(f)}\tNetCharge={total:.3f}\n")

    print(f"Wrote {out_tsv}")

if __name__ == "__main__":
    main()
