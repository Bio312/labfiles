# Error with Plot Trees and Domains R script
Replace `plotTreeAndDomains2.r` with `plotTreeAndDomains3.r`
To download:
```bash
cd ~/lab06-$MYGIT/
wget https://raw.githubusercontent.com/Bio312/labfiles/refs/heads/main/plotTreeAndDomains3.r
cd ~/lab06-$MYGIT/myoglobin
```
To run, usage is the same as for plotTreeAndDomains3.r:
```bash
Rscript --vanilla ~/lab06-$MYGIT/plotTreeAndDomains3.r <treefile> <rps-blast.out> <homologs.fas> <output.pdf>
```

# Option 1 
If you have messed up your jupyter notebook, and want the original copy of the jupyter notebook, you can grab it.
First, rename your current notebook:
```bash
mv ~/lab06-$MYGIT/myoglobin/BIO312_Myoglobin_Visualization.ipynb ~/lab06-$MYGIT/myoglobin/BIO312_Myoglobin_Visualization.old.ipynb
```

Now, grab the original version that you cloned with your lab:
```bash
cd ~/lab06-$MYGIT/myoglobin/
wget
```

# Option 2 Fixes
## Compute net charge for each PQR
For Option 2, under: *Compute net charge for each PQR*
there is a typo in the following:
```bash
out=~/lab06-$MYGIT/myoglobin/net_charges.tsv
: > "$out"   # truncate/create

for f in ~/lab06-$MYGIT/myoglobin/*.pqr; do
  base=$(basename "$f")
  Z=$(awk '$1=="ATOM"||$1=="HETATM"{s+=$9} END{printf "%.3f", s}' "$f")
  printf "%s\tNetCharge=%s\n" "$base" "$Z" >> "$out"
done
```
- The typo is that the correct column for charges is $10, not $9!
- I apologize; the output is not usable.

Here is an adjusted script that is easier to run:
Download the script:
```bash
cd ~/lab06-$MYGIT/
wget https://raw.githubusercontent.com/Bio312/labfiles/refs/heads/main/sum_pqr_charges.py
cd ~/lab06-$MYGIT/myoglobin
```
Here is how to run the script:
```bash
mamba activate bio312
python3 ~/lab06-$MYGIT/sum_pqr_charges.py ~/lab06-$MYGIT/myoglobin ~/lab06-$MYGIT/myoglobin/net_charges.tsv
```

## Fix: Compare aquatic vs terrestrial proteins
Under the section: *Compare aquatic vs terrestrial proteins*
If you have trouble running the code using cut and paste, use this python script instead. 
To get the python script:
```bash
cd ~/lab06-$MYGIT/
wget https://raw.githubusercontent.com/Bio312/labfiles/refs/heads/main/label_and_plot_net_charges.py
cd ~/lab06-$MYGIT/myoglobin
```

To run the python script:
```bash
mamba activate bio312
python3 ~/lab06-$MYGIT/label_and_plot_net_charges.py \
  --charges     ~/lab06-$MYGIT/myoglobin/net_charges.tsv \
  --species-key ~/lab06-$MYGIT/species_key.csv \
  --refseq-map  ~/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
  --out-tsv     ~/lab06-$MYGIT/myoglobin/net_charges_labeled.tsv \
  --out-plot    ~/lab06-$MYGIT/myoglobin/netcharge_boxplot.png
```

# Complete updated lab, starting at Option 2
## Option 2: Using the BASH shell

In this part of the lab, we want to **visualize myoglobin structures** and compare **electrostatic charge distributions** between aquatic and terrestrial mammals.  
Previously, we generated `.pqr` files with **pdb2pqr**. Now we’ll use simple bash + python steps to make **static PNGs** of the protein surfaces.

## Compute net charge for each PQR
Charges are already embedded in your `.pqr` files (last two columns). This script will sum them:

```bash
mamba activate bio312
python3 ~/lab06-$MYGIT/sum_pqr_charges.py ~/lab06-$MYGIT/myoglobin ~/lab06-$MYGIT/myoglobin/net_charges.tsv
```

Look at the resulting tab-separated file `net_charges.tsv` with each protein and its net charge.  
```bash
less ~/lab06-$MYGIT/myoglobin/net_charges.tsv
```

## 4. Visualize structures with charges (static PNGs)
We can’t use interactive viewers outside Jupyter, but we can generate PNG snapshots with `py3Dmol` (Python).  
First, install pymol (one time only)
```bash
mamba activate bio312
mamba install -y -c conda-forge pymol-open-source
```

First, we need to create pdbs that have our pqr charges on them:
```bash
python ~/lab06-$MYGIT/pqr_to_charges_pdb.py ~/lab06-$MYGIT/myoglobin
```
If this worked it should say "Done: created #/# charges PDBs."

Next, run the following loop to render the PDB structures colored by these charges:

```bash
mamba activate bio312
python3 ~/lab06-$MYGIT/render_charge_pngs.py ~/lab06-"$MYGIT"/myoglobin
```

**Output:**  
A `.png` image for each `.pqr` file in the same directory. These images show protein surfaces colored by charge:  
- **Blue = positive charge**  
- **Red = negative charge**  
- **White/grey = neutral**  

---

Take a look at the images. Need help figuring out which images belong to which species, and which are aquatic and terrestrial? This will build some helper files for you:
```bash
# Go to your images directory
cd ~/lab06-"$MYGIT"/myoglobin

# 1) Build helper maps:
#    - abbreviation -> status (aquatic/terrestrial)
awk -F, 'NR>1 {print $2"\t"$5}' ~/lab06-"$MYGIT"/species_key.csv \
  | sort -u > abbr_status.tsv

#    - refseq -> abbreviation (from lab03 list like Hsap|NP_...|MB)
awk -F'|' '{print $2"\t"$1}' ~/lab03-"$MYGIT"/myoglobin/myoglobin.blastp.detail.filtered.out \
  | sort -u > refseq_abbr.tsv

# 2) Map refseq -> png filename (derives refseq from the part before "__")
ls *.png 2>/dev/null \
  | sed 's/\.png$//' \
  | awk -F'__' '{print $1"\t"$0".png"}' \
  | sort -u > refseq_png.tsv

# 3) Join to get refseq -> abbr -> status
#    (join by abbreviation, so flip refseq_abbr first)
awk -F'\t' '{print $2"\t"$1}' refseq_abbr.tsv | sort -u > abbr_refseq.tsv
join -t $'\t' -1 1 -2 1 abbr_refseq.tsv abbr_status.tsv > abbr_refseq_status.tsv
awk -F'\t' '{print $2"\t"$1"\t"$3}' abbr_refseq_status.tsv | sort -u > refseq_abbr_status.tsv
# Columns now: RefSeq \t Abbr \t Status

# 4) Attach PNG filename (join on RefSeq)
join -t $'\t' -1 1 -2 1 refseq_png.tsv refseq_abbr_status.tsv > png_refseq_abbr_status.tsv
# Columns: RefSeq \t PNG \t Abbr \t Status

# 5) Split into lists and make a pretty guide
awk -F'\t' '$4=="aquatic"{print $2}' png_refseq_abbr_status.tsv > aquatic_pngs.txt
awk -F'\t' '$4=="terrestrial"{print $2}' png_refseq_abbr_status.tsv > terrestrial_pngs.txt

# Human-readable table
echo
echo "PNG GUIDE (RefSeq | PNG | Abbr | Status):"
column -t -s $'\t' png_refseq_abbr_status.tsv | tee PNG_GUIDE.txt

echo
echo "Counts -> Aquatic: $(wc -l < aquatic_pngs.txt 2>/dev/null || echo 0), Terrestrial: $(wc -l < terrestrial_pngs.txt 2>/dev/null || echo 0)"
```

Now, take a look at PNG_GUIDE.txt to figure out which PNGs are terrestrial and which are aquatic.
```
less PNG_GUIDE.txt
```
* What to pay attention to in the images
Color meaning (from .pqr charges):
- Blue: positive (Lys/Arg/His side chains are often contributors)
- Red: negative (Asp/Glu)
- White/grey: near-neutral
* Overall trend you’re testing:
Diving/semi-aquatic mammals should show more positive (blue) surface patches than terrestrial mammals. This is tied to aggregation avoidance at high intracellular myoglobin concentration.

Where the color shows up:
- Focus on the outer surface rather than the cartoon backbone (cartoon color is just a rainbow by residue index unless you changed it).
- Look for broad blue patches vs mixed/neutral surface.

Caveats:
- Our coloring reflects total per-atom charge assigned by pdb2pqr at the chosen pH; it’s a useful proxy for surface charge distribution, but not the exact “net surface charge” metric in Mirceta et al.
- Structures are not guaranteed to be aligned; orientation may differ between species.



###  Compare aquatic vs terrestrial proteins
Now combine your net charge results with the `species_key.csv` file to make a simple boxplot:

```bash
mamba activate bio312
python3 ~/lab06-$MYGIT/label_and_plot_net_charges.py \
  --charges     ~/lab06-$MYGIT/myoglobin/net_charges.tsv \
  --species-key ~/lab06-$MYGIT/species_key.csv \
  --refseq-map  ~/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
  --out-tsv     ~/lab06-$MYGIT/myoglobin/net_charges_labeled.tsv \
  --out-plot    ~/lab06-$MYGIT/myoglobin/netcharge_boxplot.png
```

**Output:**  
- A PNG file `netcharge_boxplot.png` showing net charge distributions for aquatic vs terrestrial species.  

Now, take a look at this png.


#### What you should find
- **Terrestrial species** often have lower net charge (closer to neutral).  
- **Aquatic/diving species** typically have more positive net charges.  
- The PNG images give a **qualitative** view of surface charge patches.  
- The boxplot provides a **quantitative comparison** between habitats.  

Is this what you see?

## After you have finished Option 1 or Option 2
Once you have finished the Jupyter Notebook, you should have:
- net_charge_boxplot.png
- net_charge_summary.csv

><img src="img/github.png" alt= “” width="20" height="20">[Github] These two files should be part of your GitHub repository.

Now, continue with using DSSP to extract secondary structure characteristics of myoglobin.

First, install some more software:
```
mamba activate bio312
mamba install -y -c conda-forge mdtraj pandas matplotlib
```

Josh made this script to complete the analysis of secondary structure:
```bash
mamba activate bio312
python3 ~/lab06-$MYGIT/dssp_batch_summary_mdtraj.py \
  --pdb-dir   ~/lab06-$MYGIT/myoglobin \
  --species-key ~/lab06-$MYGIT/species_key.csv \
  --refseq-map  ~/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
  --out-csv   ~/lab06-$MYGIT/myoglobin/dssp_summary.csv \
  --plots
```

This previous python script:
- Generated .dssp files
- Parse the generated .dssp files
- Counted fraction helix (H/G/I), sheet (E/B), coil (everything else)
- Reported average backbone ASA (exposedness proxy)
- Saved a CSV and make simple plots grouped by aquatic/terrestrial using species_key.csv

Examine each of these.

