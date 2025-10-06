#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd

# Use headless backend for servers/CLI
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    ap = argparse.ArgumentParser(
        description="Label PQR net charges with species info and plot by habitat."
    )
    ap.add_argument("--charges", required=True,
                    help="Input TSV with lines like: <File>\\tNetCharge=<value>")
    ap.add_argument("--species-key", required=True,
                    help="CSV with columns including 'abbreviation' and 'aquatic'")
    ap.add_argument("--refseq-map", required=True,
                    help=("Mapping file (e.g., lab03 list) where each line's first column "
                          "contains 'Abbr|RefSeq|...' (we only need Abbr and RefSeq)."))
    ap.add_argument("--out-tsv", required=True,
                    help="Output labeled TSV (adds RefSeq, Abbr, Status)")
    ap.add_argument("--out-plot", default=None,
                    help="Output PNG path for boxplot (default: alongside out-tsv)")
    return ap.parse_args()


def read_charges_table(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                fname, ztok = line.split("\t", 1)
            except ValueError:
                # malformed line; skip
                continue
            if not ztok.startswith("NetCharge="):
                continue
            try:
                z = float(ztok.split("=", 1)[1])
            except ValueError:
                continue
            rows.append((fname, z))
    if not rows:
        raise RuntimeError(f"No valid rows parsed from {path}")
    df = pd.DataFrame(rows, columns=["File", "NetCharge"])
    return df


def build_refseq_to_abbr(map_path):
    """
    Build RefSeq -> Abbr from a file where the first column (tab-separated)
    looks like 'Abbr|RefSeq|Gene'.
    """
    ref2abbr = {}
    with open(map_path) as fh:
        for raw in fh:
            raw = raw.strip()
            if not raw or raw.startswith("#"):
                continue
            first_col = raw.split("\t", 1)[0]
            parts = first_col.split("|")
            if len(parts) >= 2:
                abbr, refseq = parts[0], parts[1]
                # keep first-seen mapping
                ref2abbr.setdefault(refseq, abbr)
    if not ref2abbr:
        raise RuntimeError(f"No RefSeq->Abbr mappings parsed from {map_path}")
    return ref2abbr


def main():
    args = parse_args()

    charges_path   = os.path.expanduser(os.path.expandvars(args.charges))
    species_key    = os.path.expanduser(os.path.expandvars(args.species_key))
    refseq_map     = os.path.expanduser(os.path.expandvars(args.refseq_map))
    out_tsv        = os.path.expanduser(os.path.expandvars(args.out_tsv))
    out_plot       = (os.path.expanduser(os.path.expandvars(args.out_plot))
                      if args.out_plot else os.path.splitext(out_tsv)[0] + ".png")

    # 1) read net_charges.tsv (File \t NetCharge=VAL)
    charges = read_charges_table(charges_path)

    # 2) RefSeq accession from filename (before the "__")
    charges["RefSeq"] = charges["File"].str.split("__", n=1, expand=True)[0]

    # 3) RefSeq -> Abbr map (from lab03 list)
    ref2abbr = build_refseq_to_abbr(refseq_map)
    charges["Abbr"] = charges["RefSeq"].map(ref2abbr)

    # 4) Abbr -> Status (aquatic/terrestrial) from species_key.csv
    sk = pd.read_csv(species_key)  # expects 'abbreviation' and 'aquatic'
    if not {"abbreviation", "aquatic"}.issubset(sk.columns):
        raise RuntimeError(f"{species_key} must contain columns: abbreviation,aquatic")
    abbr2status = dict(zip(sk["abbreviation"], sk["aquatic"]))
    charges["Status"] = charges["Abbr"].map(abbr2status).fillna("unknown")

    # 5) Save labeled table
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    charges.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote {out_tsv}")

    # 6) Plot boxplot (aquatic vs terrestrial)
    df = charges[charges["Status"].isin(["aquatic", "terrestrial"])].copy()
    if df.empty:
        print("WARN: no aquatic/terrestrial rows found after mapping; "
              "check species_key.csv and the refseq mapping file.")
    else:
        ax = df.boxplot(column="NetCharge", by="Status", grid=False)
        plt.title("Net charge by habitat")
        plt.suptitle("")
        plt.ylabel("Net charge (from .pqr at your chosen pH)")
        plt.tight_layout()
        plt.savefig(out_plot, dpi=150)
        print(f"Wrote {out_plot}")

    # quick counts
    print("\nCounts by Status:")
    print(charges["Status"].value_counts(dropna=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)
