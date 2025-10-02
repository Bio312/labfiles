#!/usr/bin/env python3
import os, sys, re, glob, argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

HELIX = set(['H','G','I'])
SHEET = set(['E','B'])

def log(msg, stream=sys.stderr):
    print(msg, file=stream, flush=True)

def ensure_dir(p):
    if p:
        os.makedirs(p, exist_ok=True)

def load_refseq_to_abbr(map_path: str):
    """Parse lines like: Abbr|RefSeq|Gene (in the first column)."""
    ref2abbr = {}
    with open(map_path) as fh:
        for line in fh:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            first = line.split('\t',1)[0]
            parts = first.split('|')
            if len(parts) >= 2:
                abbr, refseq = parts[0], parts[1]
                ref2abbr.setdefault(refseq, abbr)
    return ref2abbr

def prefer_alphaFold(pdb_paths):
    for p in pdb_paths:
        if '__AF-' in os.path.basename(p):
            return p
    return pdb_paths[0]

def summarize_one(pdb_path: str):
    """
    Use mdtraj to:
      - load structure
      - compute DSSP per residue
      - compute SASA per residue
    Return dict with counts/fractions + mean_ASA.
    """
    # mdtraj wants a topology & coordinates; PDB is fine
    traj = md.load(pdb_path)

    # compute_dssp returns array of DSSP codes (per residue)
    dssp = md.compute_dssp(traj, simplified=False)[0]  # first (and only) frame
    # Map to helix/sheet/coil buckets
    n_total = len(dssp)
    n_helix = sum(1 for s in dssp if s in HELIX)
    n_sheet = sum(1 for s in dssp if s in SHEET)
    n_coil  = n_total - n_helix - n_sheet

    # SASA: per-atom → per-residue mean
    # mdtraj.shrake_rupley returns per-atom SASA (n_atoms) for each frame
    sasa_atom = md.shrake_rupley(traj)  # shape: (n_frames, n_atoms)
    sasa_atom = sasa_atom[0]            # 1 frame
    # Aggregate to residues
    top = traj.topology
    res_sasa = []
    for res in top.residues:
        atom_idx = [a.index for a in res.atoms]
        if atom_idx:
            res_sasa.append(float(np.mean(sasa_atom[atom_idx])))
    mean_ASA = float(np.mean(res_sasa)) if res_sasa else 0.0

    return dict(
        n_total=n_total,
        n_helix=n_helix,
        n_sheet=n_sheet,
        n_coil=n_coil,
        frac_helix=n_helix/n_total if n_total else 0.0,
        frac_sheet=n_sheet/n_total if n_total else 0.0,
        frac_coil=n_coil/n_total if n_total else 0.0,
        mean_ASA=mean_ASA
    )

def make_plots(df: pd.DataFrame, out_dir: str, by_col: str='status'):
    ensure_dir(out_dir)
    for col in ['frac_helix','frac_sheet','frac_coil','mean_ASA']:
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
    ap = argparse.ArgumentParser(description="Summarize helix/sheet/coil + ASA via mdtraj (no mkdssp).")
    ap.add_argument('--pdb-dir', required=True, help='Directory with input .pdb files')
    ap.add_argument('--species-key', required=True, help='CSV with "abbreviation" and "aquatic"')
    ap.add_argument('--refseq-map', required=True, help='File mapping RefSeq to species abbreviation')
    ap.add_argument('--out-csv', required=True, help='Summary CSV path')
    ap.add_argument('--plots', action='store_true', help='Also write boxplots next to CSV')
    args = ap.parse_args()

    pdb_dir  = os.path.expanduser(os.path.expandvars(args.pdb_dir))
    out_csv  = os.path.expanduser(os.path.expandvars(args.out_csv))
    key_csv  = os.path.expanduser(os.path.expandvars(args.species_key))
    map_path = os.path.expanduser(os.path.expandvars(args.refseq_map))

    # Load species key
    key = pd.read_csv(key_csv)
    abbr2status = dict(zip(key['abbreviation'], key['aquatic']))

    # Load RefSeq->abbr map
    ref2abbr = load_refseq_to_abbr(map_path)
    if not ref2abbr:
        log(f"WARNING: no mappings parsed from {map_path}")

    # Collect PDBs and group by RefSeq (prefix before '__')
    pdb_files = sorted(glob.glob(os.path.join(pdb_dir, '*.pdb')))
    if not pdb_files:
        log(f"No .pdb files found in {pdb_dir}")
        sys.exit(1)

    refseq_to_pdbs = {}
    for p in pdb_files:
        refseq = os.path.basename(p).split('__')[0]
        refseq_to_pdbs.setdefault(refseq, []).append(p)

    records = []
    for refseq, paths in sorted(refseq_to_pdbs.items()):
        chosen = prefer_alphaFold(paths)
        base = os.path.basename(chosen)
        try:
            metrics = summarize_one(chosen)
        except Exception as e:
            log(f"Skipping {base}: {e}")
            continue

        abbr   = ref2abbr.get(refseq, 'unknown')
        status = abbr2status.get(abbr, 'unknown')

        metrics['RefSeq'] = refseq
        metrics['abbr']   = abbr
        metrics['status'] = status
        metrics['source'] = 'AF' if '__AF-' in base else 'SWM'
        records.append(metrics)

    if not records:
        log("No summaries produced.")
        sys.exit(3)

    df = pd.DataFrame.from_records(records)
    ensure_dir(os.path.dirname(out_csv) or '.')
    df.to_csv(out_csv, index=False)
    print(f"Wrote: {out_csv}")
    print(df.head())

    if args.plots:
        plot_dir = os.path.dirname(out_csv) or '.'
        make_plots(df, plot_dir, by_col='status')

if __name__ == '__main__':
    main()
