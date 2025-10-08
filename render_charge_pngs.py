#!/usr/bin/env python3
import os
import sys
import glob

# Use PyMOL's embedded, headless API
try:
    import pymol2
except ImportError as e:
    sys.stderr.write("ERROR: Could not import pymol2. Install with:\n"
                     "  mamba install -y -c conda-forge pymol-open-source\n")
    sys.exit(1)

def render_one(pdb_path, width=1200, height=900, dpi=200, zoom=5):
    """
    Load a .charges.pdb (B-factor holds per-atom charge), show surface,
    color by charge (B-factor), and save a PNG next to the file.
    """
    base, _ = os.path.splitext(pdb_path)
    out_png = base + ".png"

    # Start a fresh headless PyMOL instance for each file (keeps things clean)
    with pymol2.PyMOL() as pm:
        cmd = pm.cmd
        pm.start()  # start the engine

        # Load and set up scene
        cmd.load(pdb_path, "m")
        cmd.hide("everything", "m")
        cmd.show("surface", "m")

        # Color the surface by B-factor (which stores atomic charge from the PQR)
        # red=negative, white=0, blue=positive
        cmd.spectrum("b", "red_white_blue", "m")

        cmd.bg_color("white")
        cmd.orient("m")
        cmd.zoom("m", zoom)
        cmd.set("ray_opaque_background", 0)

        # Render a nice PNG; 'ray=1' forces raytracing for quality
        # NOTE: in API, pass width/height to png; dpi is stored in metadata
        cmd.png(out_png, width=width, height=height, dpi=dpi, ray=1)

    # Verify the PNG was created
    if not os.path.exists(out_png) or os.path.getsize(out_png) == 0:
        raise RuntimeError(f"PNG not created: {out_png}")
    return out_png

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python3 render_charge_pngs.py <directory_with_charges_pdb>\n")
        sys.exit(2)

    in_dir = os.path.abspath(os.path.expanduser(os.path.expandvars(sys.argv[1])))
    if not os.path.isdir(in_dir):
        sys.stderr.write(f"ERROR: Not a directory: {in_dir}\n")
        sys.exit(2)

    pattern = os.path.join(in_dir, "*.charges.pdb")
    files = sorted(glob.glob(pattern))
    if not files:
        sys.stderr.write(f"No *.charges.pdb files found in {in_dir}\n")
        sys.exit(1)

    print(f"Found {len(files)} charges PDBs")
    ok = 0
    for f in files:
        rel = os.path.basename(f)
        try:
            print(f"Rendering {rel} ...", end="", flush=True)
            out = render_one(f)
            print(f" OK -> {os.path.basename(out)}")
            ok += 1
        except Exception as e:
            print(f" FAIL\n  Reason: {e}")

    print(f"\nDone. {ok}/{len(files)} images created.")

if __name__ == "__main__":
    main()
