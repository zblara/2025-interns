#!/bin/bash

SECONDS=0

echo -e "#####################################################################################

  ZBL Ｍｅｔａｇｅｎｏｍｉｃｓ Ｐｒｅ－ｐｒｏｃｅｓｓｉｎｇ Ｐｉｐｅｌｉｎｅ

  This bash script is optimized for concatenating paired-end filtered reads.

  Usage:
	Enter PATH to directory where raw reads are stored.
	Indicate one folder up. 
	
  Example:
	PE reads: /media/zbl/Storage/mgs(/raw)
	
  Output:
  	1. Concatenated reads in ~/path/filtered_concatenated.

  ©Z.B. Lara, 2024

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

# No output to this folder. For further grouping by user only.

# Function to concatenate paired-end reads
concatenate_reads() {
    local R1_file="$1"
    local DIR="$2"

    # Extract BASE from the R1 filename
    BASE=$(basename "$R1_file" | cut -d '_' -f 1,2)

    # Construct the corresponding R2 filename
    R2_file="${R1_file/_R1_/_R2_}"

    # Check if the corresponding R2 file exists
    if [ -e "$R2_file" ]; then
    
        echo -e "\nConcatenating paired end reads for $BASE.\n" $(date -u)
        
        # Define the output filenames for mapping
        filtered_1="${DIR}/filter_human/${BASE}_filtered_R1.fastq.gz"
        filtered_2="${DIR}/filter_human/${BASE}_filtered_R2.fastq.gz"
        concatenated="${DIR}/filter_human/${BASE}_filtered_concatenated.fastq.gz"
        
    		check_file() {
    		if [ -f "$1" ]; then
        		echo "Skipping $2: File $1 already exists."
        		return 1
    		else
        		echo "File $1 does not exist."
        		return 0
    		fi
			}
        
        # Check if mapping i/o files already exist
        if ! check_file "$concatenated" "concatenated output $concatenated" ; then
            continue
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

# Find all R1 files and pass them to xargs to run in parallel
find "${DIR}/filter_human/" -type f -name "*_R1.fastq.gz" | xargs -I {} -P "4" bash -c 'concatenate_reads "$@"' _ "{}" "${DIR}"

# After the script finishes, use "cat" if you wish to further group files.
# Example: "cat A1_1_filtered_concatenated.fastq.gz A2_1_filtered_concatenated.fastq.gz .... \
#			> normW0_filtered_concatenated.fastq.gz"
mkdir -p "${DIR}/filtered_concatenated/"

echo -e "\nPipeline finished." $(date -u)
echo "Elapsed time: $SECONDS seconds"
echo -e "Thank you for your patience.\n"
