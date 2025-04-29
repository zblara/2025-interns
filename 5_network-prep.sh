!/bin/bash

SECONDS=0

# ZBL Bipartite Network Preparation Script
# Purpose: Build a taxon-AMR bipartite edge list from Bracken and AMRFinderPlus outputs
# Date: 2025-04-29
# Author: Z.B. Lara

# ---------------------------
# USER CONFIGURATION SECTION
# ---------------------------
BRACKEN_DIR="./bracken_outputs"           # Directory containing Bracken report files
AMR_DIR="./amrfinder_outputs"             # Directory containing AMRFinderPlus .tsv outputs
OUTDIR="./amr_annotation/bipartite_network" # Output directory
EDGE_LIST="$OUTDIR/edge_list.tsv"         # Final edge list file

mkdir -p "$OUTDIR"
echo -e "Sample\tTaxon\tAMR_Gene" > "$EDGE_LIST"

# ---------------------------
# MAIN PROCESSING LOOP
# ---------------------------

for bracken_file in "$BRACKEN_DIR"/*.tsv; do
    sample=$(basename "$bracken_file" .tsv)
    amr_file="$AMR_DIR/${sample}_amr.tsv"

    if [[ ! -f "$amr_file" ]]; then
        echo "Warning: AMR file not found for $sample. Skipping." >&2
        continue
    fi

    # Extract top taxa from Bracken (e.g., species-level)
    awk -F'\t' 'NR>1 && $3=="S" && $5 > 0 {print $1}' "$bracken_file" > "$OUTDIR/${sample}_taxa.tmp"

    # Extract AMR genes
    awk -F'\t' 'NR>1 && $4=="AMR" {print $2}' "$amr_file" | sort -u > "$OUTDIR/${sample}_amr.tmp"

    # Join every taxon with every AMR gene (co-occurrence assumption)
    while read taxon; do
        while read gene; do
            echo -e "$sample\t$taxon\t$gene" >> "$EDGE_LIST"
        done < "$OUTDIR/${sample}_amr.tmp"
    done < "$OUTDIR/${sample}_taxa.tmp"

done

# Clean up temporary files
rm "$OUTDIR"/*.tmp

echo "Bipartite edge list written to $EDGE_LIST"

echo -e "\n Network preparation pipeline finished." $(date -u)
echo "Elapsed time: $SECONDS seconds"
echo -e "Thank you for your patience.\n"
