#!/bin/bash

SECONDS=0

echo -e "#####################################################################################

  ZBL Metagenomics Pre-processing Pipeline

  This bash script is optimized for concatenating paired-end filtered reads.

  Usage:
	Enter PATH to directory where raw reads are stored.
	Indicate one folder up. 
	
  Example:
	PE reads: /media/zbl/Storage/mgs(/raw)

  Output:
  	1. Concatenated reads in ~/path/filtered_concatenated.

  Date: 2025-04-29
  Author: Z.B. Lara

######################################################################################"

# Function to check file existence and proceed accordingly
check_file() {
    if [ -f "$1" ]; then
        echo "Skipping $2: File $1 already exists."
        return 1
    else
        echo "File $1 does not exist."
        return 0
    fi
}

# Prompt for the directory path
read -p "Enter PATH to files: " DIR

# Check if the directory exists
if [ ! -d "$DIR" ]; then
    echo -e "Error: The directory does not exist.\n"
    exit 1
fi

# Create output directories if they don't exist
mkdir -p "${DIR}/filtered_concatenated/"
mkdir -p "${DIR}/filter_human/"

# Function to concatenate paired-end reads
concatenate_reads() {
    local R1_file="$1"
    local DIR="$2"

    # Extract BASE from the R1 filename
    BASE=$(basename "$R1_file" | cut -d '_' -f 1)

    # Construct the corresponding R2 filename
    R2_file="${R1_file/_R1_/_R2_}"

    # Check if the corresponding R2 file exists
    if [ -e "$R2_file" ]; then
        
        echo -e "\nConcatenating paired end reads for $BASE.\n" $(date -u)
        
        # Define the output filenames for concatenation
        filtered_1="${DIR}/filter_human/${BASE}_filtered_R1.fastq.gz"
        filtered_2="${DIR}/filter_human/${BASE}_filtered_R2.fastq.gz"
        concatenated="${DIR}/filtered_concatenated/${BASE}_filtered_concatenated.fastq.gz"
        
        # Check if concatenated output file already exists
        if ! check_file "$concatenated" "concatenated output $concatenated"; then
            return  # Exit early if the file already exists
        fi
        
        # Execute the command for concatenating paired-end reads
        paste <(zcat "${filtered_1}") <(zcat "${filtered_2}") \
        | awk 'NR%4==1{gsub("_filtered", "_concatenated", $1)} 1' \
        | gzip > "${concatenated}"
        
        # Add verbose option to track progress
        echo "Processed concatenating of $BASE." $(date -u)
    else
        echo "Corresponding R2 file not found for $R1_file."
    fi
}

export -f concatenate_reads
export -f check_file

# Find all R1 files and pass them to xargs to run in parallel
find "${DIR}/filter_human/" -type f -name "*_R1.fastq.gz" | xargs -I {} -P 4 bash -c 'concatenate_reads "$@"' _ "{}" "${DIR}"

echo -e "\nPipeline finished." $(date -u)
echo "Elapsed time: $SECONDS seconds"
echo -e "Thank you for your patience.\n"
