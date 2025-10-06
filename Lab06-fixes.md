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
wget https://raw.githubusercontent.com/Bio312/labfiles/refs/heads/main/sum_pqr_charges.py
cd ~/lab06-$MYGIT/
wget
cd cd ~/lab06-$MYGIT/myoglobin
```
Here is how to run the script:
```bash
python3 sum_pqr_charges.py ~/lab06-$MYGIT/myoglobin ~/lab06-$MYGIT/myoglobin/net_charges.tsv
```

## Fix: Compare aquatic vs terrestrial proteins
Under the section: *Compare aquatic vs terrestrial proteins*
If you have trouble running the code using cut and paste, use this python script instead. 
To get the python script:
```wget
cd ~/lab06-$MYGIT/
wget
cd cd ~/lab06-$MYGIT/myoglobin
```

To run the python script:
```R
python3 label_and_plot_net_charges.py \
  --charges     ~/lab06-$MYGIT/myoglobin/net_charges.tsv \
  --species-key ~/lab06-$MYGIT/species_key.csv \
  --refseq-map  ~/lab03-$MYGIT/myoglobin/myoglobin.blastp.detail.filtered.out \
  --out-tsv     ~/lab06-$MYGIT/myoglobin/net_charges_labeled.tsv \
  --out-plot    ~/lab06-$MYGIT/myoglobin/netcharge_boxplot.png
```

