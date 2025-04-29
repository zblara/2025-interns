#!/bin/bash

SECONDS=0

echo -e "#####################################################################################

  ZBL AMR Annotation Pipeline

  This bash script is optimized for AMR gene annotation using AMRFinderPlus
  on metagenomic reads in FASTQ format.

  Input: Concatenated PE reads (e.g., *_filtered_concatenated.fastq.gz)

  Output: 
    1. AMRFinderPlus results in ~/path/amr_annotation/

  © Z.B. Lara, 2024

######################################################################################"

# Function to check file existence
check_file() {
    if [ -f "$1" ]; then
        echo "Skipping $2: File $1 already exists."
        return 1
    else
        echo "File $1 does not exist. Proceeding."
        return 0
    fi
}

# Function to check if file was successfully created
check_success() {
    if [ -f "$1" ]; then
        echo -e "\nFile $1 processing successful." $(date -u)
        return 1
    else
        echo -e "❌ Error: Processing $1 not successful.\n"
        return 0
    fi
}

# Prompt for the directory path
read -p "Enter PATH to directory with filtered_concatenated reads: " DIR

# Check if the directory exists
if [ ! -d "$DIR" ]; then
    echo -e "❌ Error: The directory does not exist.\n"
    exit 1
fi

# Create output directory
mkdir -p "${DIR}/amr_annotation/"

# Activate conda environment
echo -e "\nActivating metagenomics conda environment."
source ~/miniconda3/bin/activate mgs
conda env list

# Loop through all concatenated filtered files
for fq_file in "${DIR}/filtered_concatenated/"*_filtered_concatenated.fastq.gz; do
    BASE=$(basename "$fq_file" _filtered_concatenated.fastq.gz)

    echo -e "\n[INFO] Starting AMR annotation for $BASE. $(date -u)"

    output_file="${DIR}/amr_annotation/${BASE}.amrfinder.txt"

    # Skip if file already exists
    if ! check_file "$output_file" "AMRFinderPlus output for $BASE"; then
        continue
    fi

    # Run AMRFinderPlus
    amrfinder -n "$fq_file" \
              -o "$output_file" \
              --organism other \
              --threads 8 \
              --plus \
              --name "${BASE}"

    # Check for successful output
    if ! check_success "$output_file"; then
        continue
    fi

    echo "Finished AMRFinderPlus for $BASE. $(date -u)"
done

echo -e "\n Pipeline finished. $(date -u)"
echo "Elapsed time: $SECONDS seconds"
echo -e "Thank you for your patience.\n"

