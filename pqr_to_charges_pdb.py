#!/usr/bin/env python3
import os, sys, glob, re

HELP = """\
Usage: pqr_to_charges_pdb.py <input_dir_with_pqr>
Creates *.charges.pdb next to each .pqr by copying coords and putting per-atom charge into the B-factor column.
"""

# Regex: capture the last two floats on a line (charge and radius in PQR)
TAIL_FLOATS = re.compile(r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$')

def convert_pqr_to_charges_pdb(pqr_path, out_path):
    """
    Read PQR, write PDB with charge stored in B-factor column.
    PDB columns (fixed width):
      - x,y,z: cols 31-54
      - occupancy: cols 55-60
      - b-factor: cols 61-66  <-- we put charge here
    """
    out = []
    with open(pqr_path) as fh:
        for line in fh:
            if line.startswith(('ATOM', 'HETATM')):
                m = TAIL_FLOATS.search(line)
                if not m:
                    continue
                charge = float(m.group(1))  # charge is the first of the last two floats
                # Expand PQR fields to PDB-like fixed columns
                # We’ll rewrite occ to 1.00 and b-factor to charge
                rec   = line[:6]                      # ATOM/HETATM (cols 1-6)
                rest  = line[6:66]                    # up through b-factor area
                # Re-extract XYZ safely from PQR to preserve positions
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    # Fallback: token parse (some PQRs don’t align perfectly)
                    toks = line.split()
                    # last two are charge,radius; xyz are toks[-5:-2]
                    x, y, z = map(float, toks[-5:-2])

                # Build a fresh PDB ATOM/HETATM line with standard formatting
                # We keep atom names/residue/chain/resi from the original slices
                # and inject new XYZ, occ, b-factor
                new = (
                    f"{line[:6]}"                               # ATOM/HETATM
                    f"{line[6:30]}"                             # atom serial/name/alt/ res/chain/resi/icode
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"                 # xyz
                    f"{1.00:6.2f}{charge:6.2f}"                 # occ=1.00, b=charge
                    f"{line[66:].rstrip()}\n"                   # element/charge if present
                )
                out.append(new)
            elif line.startswith(('TER','END','ENDMDL','MODEL','REMARK','#')):
                # Copy safe meta lines through
                out.append(line if line.endswith('\n') else line+'\n')

    # Ensure file ends with END
    if not any(l.startswith('END') for l in out[-3:]):
        out.append("END\n")

    with open(out_path, 'w') as w:
        w.writelines(out)

def main():
    if len(sys.argv) != 2:
        print(HELP, file=sys.stderr)
        sys.exit(2)

    in_dir = os.path.abspath(os.path.expanduser(os.path.expandvars(sys.argv[1])))
    if not os.path.isdir(in_dir):
        print(f"ERROR: not a directory: {in_dir}", file=sys.stderr)
        sys.exit(2)

    pqr_files = sorted(glob.glob(os.path.join(in_dir, "*.pqr")))
    if not pqr_files:
        print(f"No .pqr files found in {in_dir}", file=sys.stderr)
        sys.exit(1)

    made = 0
    for pqr in pqr_files:
        out_pdb = os.path.splitext(pqr)[0] + ".charges.pdb"
        try:
            convert_pqr_to_charges_pdb(pqr, out_pdb)
            print(f"OK  {os.path.basename(out_pdb)}")
            made += 1
        except Exception as e:
            print(f"FAIL {os.path.basename(pqr)}  ({e})", file=sys.stderr)
    print(f"Done: created {made}/{len(pqr_files)} charges PDBs.")

if __name__ == "__main__":
    main()
