#!/usr/bin/env python3
"""
Batch DSSP summary for AlphaFold/SwissModel PDBs (with species mapping).

Fixes:
- Uses a supplied mapping file to map RefSeq -> species abbreviation
  (e.g., myoglobin.blastp.detail.filtered.out with lines like
   Abbr|RefSeq|Gene in the first column).
- Joins species abbreviation to species_key.csv to get aquatic/terrestrial.
- De-duplicates multiple PDBs per RefSeq, preferring AlphaFold over SWISS-MODEL.

Outputs:
- CSV with n_total, n_helix, n_sheet, n_coil, fractions, mean_ASA, RefSeq, abbr, status
- Optional boxplots by status.

Example:
  python3 dssp_batch_summary.py \
    --pdb-dir   ~/lab06-$MYGIT/myoglobin/AF_structs \
    --dssp-dir  ~/lab06-$MYGIT/myoglobin/AF_structs/dssp_out \
    --species-key ~/lab06-$MYGIT/species_key.csv \
    --refseq-map  ~/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
    --out-csv   ~/lab06-$MYGIT/myoglobin/AF_structs/dssp_summary.csv \
    --plots
"""

import os
import sys
import re
import glob
import argparse
import subprocess
import shutil
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HELIX = set(['H','G','I'])
SHEET = set(['E','B'])

def log(msg: str, stream=sys.stderr):
    print(msg, file=stream, flush=True)

def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)

def clean_pdb(in_pdb: str, out_pdb: str):
    """
    Keep only ATOM/HETATM/TER/END lines to avoid DSSP parser issues.
    """
    with open(in_pdb) as fin, open(out_pdb, 'w') as fout:
        for line in fin:
            if line.startswith(('ATOM', 'HETATM', 'TER', 'END')):
                fout.write(line)

def run_mkdssp(clean_pdb_path: str, out_dssp_path: str) -> bool:
    mkdssp = shutil.which("mkdssp")
    if mkdssp is None:
        log("ERROR: mkdssp not found in PATH.")
        return False
    try:
        subprocess.run([mkdssp, clean_pdb_path, out_dssp_path],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError as e:
        log(f"mkdssp failed for {clean_pdb_path} -> {out_dssp_path}")
        err = e.stderr.decode(errors="ignore")
        if err:
            log(err)
        return False

def parse_dssp_counts(dssp_path: str):
    """
    Parse mkdssp table for secondary structure and ASA.
    """
    in_table = False
    n_total = n_helix = n_sheet = n_coil = 0
    asas = []

    with open(dssp_path) as fh:
        for line in fh:
            if line.startswith('  #  RESIDUE'):
                in_table = True
                continue
            if not in_table:
                continue
            if len(line) < 50:
                continue

            ss = line[16].strip()  # secondary structure char
            asa_str = line[34:38].strip()  # ASA (classic pos kept by mkdssp)
            try:
                asa = float(asa_str) if asa_str else 0.0
            except ValueError:
                asa = 0.0

            n_total += 1
            if ss in HELIX:
                n_helix += 1
            elif ss in SHEET:
                n_sheet += 1
            else:
                n_coil += 1
            asas.append(asa)

    if n_total == 0:
        raise AssertionError(f"No residues parsed from {dssp_path}")

    return {
        'n_total': n_total,
        'n_helix': n_helix,
        'n_sheet': n_sheet,
        'n_coil':  n_coil,
        'frac_helix': n_helix / n_total,
        'frac_sheet': n_sheet / n_total,
        'frac_coil':  n_coil  / n_total,
        'mean_ASA':   (sum(asas)/len(asas)) if asas else 0.0,
    }

def load_refseq_to_abbr(map_path: str):
    """
    Build a dict RefSeq -> species abbreviation from a file whose lines look like:
      Abbr|RefSeq|Gene
    (We only need the Abbr (field1) and RefSeq (field2).)
    """
    ref2abbr = {}
    with open(map_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):  # skip blanks/comments
                continue
            # lines might be TAB-delimited columns; the first column holds "Abbr|RefSeq|Gene"
            first = line.split('\t', 1)[0]
            parts = first.split('|')
            if len(parts) >= 2:
                abbr, refseq = parts[0], parts[1]
                # take first seen mapping
                ref2abbr.setdefault(refseq, abbr)
    return ref2abbr

def prefer_alphaFold(pdb_paths):
    """
    From a list of PDB paths for the same RefSeq, pick AlphaFold if present,
    else pick the first.
    """
    for p in pdb_paths:
        if '__AF-' in os.path.basename(p):
            return p
    return pdb_paths[0]

def make_plots(df: pd.DataFrame, out_dir: str, by_col: str = 'status'):
    ensure_dir(out_dir)
    for col in ['frac_helix', 'frac_sheet', 'frac_coil', 'mean_ASA']:
        ax = df.boxplot(column=col, by=by_col, grid=False)
        ax.set_title(f'{col} by {by_col}')
        ax.set_ylabel(col)
        plt.suptitle('')
        plt.tight_layout()
        out_png = os.path.join(out_dir, f'{col}_by_{by_col}.png')
        plt.savefig(out_png, dpi=150)
        plt.close()
        print(f"Wrote plot: {out_png}")

def main():
    ap = argparse.ArgumentParser(description="Run DSSP on PDBs and summarize helix/sheet/coil + ASA (with species mapping).")
    ap.add_argument('--pdb-dir', required=True, help='Directory containing input .pdb files')
    ap.add_argument('--dssp-dir', required=True, help='Directory to write cleaned PDBs and .dssp files')
    ap.add_argument('--species-key', required=True, help='CSV with headers including "abbreviation" and "aquatic"')
    ap.add_argument('--refseq-map', required=True, help='File mapping RefSeq to species abbreviation (e.g., myoglobin.blastp.detail.filtered.out)')
    ap.add_argument('--out-csv', required=True, help='Path to summary CSV to write')
    ap.add_argument('--plots', action='store_true', help='If set, write boxplots alongside out-csv')
    args = ap.parse_args()

    pdb_dir   = os.path.expanduser(os.path.expandvars(args.pdb_dir))
    dssp_dir  = os.path.expanduser(os.path.expandvars(args.dssp_dir))
    out_csv   = os.path.expanduser(os.path.expandvars(args.out_csv))
    key_csv   = os.path.expanduser(os.path.expandvars(args.species_key))
    map_path  = os.path.expanduser(os.path.expandvars(args.refseq_map))

    clean_dir = os.path.join(dssp_dir, 'clean')
    ensure_dir(dssp_dir)
    ensure_dir(clean_dir)

    # Load species key
    try:
        key = pd.read_csv(key_csv)
        abbr2status = dict(zip(key['abbreviation'], key['aquatic']))
    except Exception as e:
        log(f"ERROR reading species key {key_csv}: {e}")
        sys.exit(2)

    # Load RefSeq -> Abbr mapping
    if not os.path.exists(map_path):
        log(f"ERROR: refseq-map file not found: {map_path}")
        sys.exit(2)
    ref2abbr = load_refseq_to_abbr(map_path)
    if not ref2abbr:
        log(f"WARNING: Could not parse any RefSeq->abbr mappings from {map_path}")

    # Group PDBs by RefSeq (prefix before __)
    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, '*.pdb')))
    if not pdb_files:
        log(f"No .pdb files found in {pdb_dir}.")
        sys.exit(1)

    refseq_to_pdbs = {}
    for p in pdb_files:
        base = os.path.basename(p)
        refseq = base.split('__')[0]
        refseq_to_pdbs.setdefault(refseq, []).append(p)

    # Warn if CIF dictionaries likely missing (only warn; you fetched earlier)
    libcifpp_dir = os.environ.get('LIBCIFPP_DATA_DIR')
    if not libcifpp_dir or not os.path.exists(os.path.join(libcifpp_dir, 'components.cif')):
        log("WARNING: LIBCIFPP_DATA_DIR not set or components.cif not found. "
            "If mkdssp fails, fetch the wwPDB dictionaries and export LIBCIFPP_DATA_DIR.\n", stream=sys.stderr)

    # Process: choose one PDB per RefSeq (prefer AF), clean, run DSSP, parse
    records = []
    for refseq, paths in sorted(refseq_to_pdbs.items()):
        chosen_pdb = prefer_alphaFold(paths)
        base = os.path.basename(chosen_pdb)

        clean_path = os.path.join(clean_dir, base.replace('.pdb', '.clean.pdb'))
        dssp_path  = os.path.join(dssp_dir, base.replace('.pdb', '.dssp'))

        try:
            clean_pdb(chosen_pdb, clean_path)
        except Exception as e:
            log(f"Skipping (clean failed): {base} :: {e}")
            continue

        if not os.path.exists(dssp_path):
            ok = run_mkdssp(clean_path, dssp_path)
            if not ok:
                log(f"Skipping (mkdssp failed): {base}")
                continue

        try:
            metrics = parse_dssp_counts(dssp_path)
        except Exception as e:
            log(f"Skipping (parse failed): {base} :: {e}")
            continue

        abbr = ref2abbr.get(refseq, 'unknown')
        status = abbr2status.get(abbr, 'unknown')

        metrics['RefSeq'] = refseq
        metrics['abbr']   = abbr
        metrics['status'] = status
        metrics['source'] = 'AF' if '__AF-' in base else 'SWM'
        records.append(metrics)

    if not records:
        log("No DSSP summaries produced. Exiting.")
        sys.exit(3)

    df = pd.DataFrame.from_records(records)
    ensure_dir(os.path.dirname(out_csv) or '.')
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")
    print(df.head())

    if args.plots:
        plot_dir = os.path.dirname(out_csv) or '.'
        make_plots(df, plot_dir, by_col='status')

if __name__ == '__main__':
    main()
