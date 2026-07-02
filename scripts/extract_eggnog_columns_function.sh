#!/usr/bin/env bash

# Usage: ./extract_function.sh COLUMN_NAME
# Example: ./extract_function.sh KEGG_ko
# Will produce: metagenomic_browning_KEGG_ko.tsv, metagenomic_browning_KEGG_ko_expanded.tsv, metagenomic_browning_KEGG_ko_counts.tsv

# Column map (adjust if eggNOG output changes):
# query = 1
# seed_ortholog = 2
# evalue = 3
# score = 4
# eggNOG_OGs = 5
# max_annot_lvl = 6
# COG_category = 7
# Description = 8
# Preferred_name = 9
# GOs = 10
# EC = 11
# KEGG_ko = 12
# KEGG_Pathway = 13
# KEGG_Module = 14
# KEGG_Reaction = 15
# KEGG_rclass = 16
# BRITE = 17
# KEGG_TC = 18
# CAZy = 19
# BiGG_Reaction = 20
# PFAMs = 21

declare -A colmap
colmap[COG_category]=7
colmap[GOs]=10
colmap[EC]=11
colmap[KEGG_ko]=12
colmap[KEGG_Pathway]=13
colmap[KEGG_Module]=14
colmap[KEGG_Reaction]=15
colmap[KEGG_rclass]=16
colmap[BRITE]=17
colmap[KEGG_TC]=18
colmap[CAZy]=19
colmap[BiGG_Reaction]=20
colmap[PFAMs]=21

colname=$1
colnum=${colmap[$colname]}

if [ -z "$colnum" ]; then
    echo "Error: Column $colname not recognized."
    echo "Available: ${!colmap[@]}"
    exit 1
fi

echo "Extracting column $colname (col $colnum)..."

# Combine all samples into one file
outfile="metagenomic_browning_${colname}.tsv"
echo -e "Sample\tQuery\t${colname}" > "$outfile"
for f in *_eggnog.emapper.annotations; do
    sample=$(basename "$f" _eggnog.emapper.annotations)
    awk -F'\t' -v s="$sample" -v c="$colnum" 'BEGIN{OFS="\t"} !/^#/ && $c != "-" {print s, $1, $c}' "$f"
done >> "$outfile"

# Expand multiple entries (if comma-separated)
expanded="metagenomic_browning_${colname}_expanded.tsv"
awk -F'\t' -v c=3 'NR==1{print; next}
  $c!="-"{gsub(/ /,"",$c); n=split($c,a,",");
  for(i=1;i<=n;i++) print $1"\t"$2"\t"a[i]}' OFS='\t' "$outfile" > "$expanded"

# Make counts table
counts="metagenomic_browning_${colname}_counts.tsv"
awk -F'\t' 'NR>1{print $1"\t"$3}' OFS='\t' "$expanded" \
  | sort \
  | uniq -c \
  | awk '{print $2"\t"$3"\t"$1}' OFS='\t' > "$counts"

echo "Done! Files created:"
echo " - $outfile"
echo " - $expanded"
echo " - $counts"
