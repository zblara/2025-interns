#!/bin/bash

SECONDS=0

echo -e "#####################################################################################

  ZBL Ｍｅｔａｇｅｎｏｍｉｃｓ Ｃｌａｓｓｉｆｉｃａｔｉｏｎ Ｐｉｐｅｌｉｎｅ

  This bash script is optimized for the taxonomic classification of filtered reads.

  Usage:
	Enter Kraken2 DATABASE name and PATH to directory where filtered reads are stored.
	Indicate one folder up. 
	
  Example:
	Enter DATABASE name (sample): standard8
	PE reads: /media/zbl/Storage/mgs(/filtered)
	
  Output:
  	1. Kraken2 output in ~/path/kraken_fc_DBname/
  	2. Bracken output in ~/path/bracken_fc_DBname/

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

check_success() {
    if [ -f "$1" ]; then
        echo -e "\nFile $1 processing successful." $(date -u)
        return 1
    else
        echo "Error: Processing $1 not successful."
        return 0
    fi
}

# Prompt for the database name
read -p "Enter Kraken DATABASE name (db): " DB

# Prompt for the directory path
read -p "Enter PATH to files: " DIR

# Check if the directory exists
if [ ! -d "$DIR" ]; then
    echo -e "Error: The directory does not exist.\n"
    exit 1
fi

# Create output directories if they don't exist
mkdir -p "${DIR}/classification/kraken2_fc_${DB}/"
mkdir -p "${DIR}/classification/bracken_fc_${DB}/"
mkdir -p "${DIR}/classification/pavian_fc_${DB}/"

# Activate conda environment
echo -e "\nActivating metagenomics conda environment."
source ~/miniconda3/bin/activate mgs
conda env list

# Loop through all R1 files matching the pattern *_R1_001.fastq.gz
for file in "${DIR}/filtered_concatenated/"*_filtered_concatenated.fastq.gz; do
    # Extract BASE from the R1 filename
    BASE=$(basename "$file" | cut -d '_' -f 1)
    
        echo -e "\nAttempting to perform alignment.\n" $(date -u)
        
        # Define the output filenames for mapping
		kraken2_output="${DIR}/classification/kraken2_fc_${DB}/${BASE}.kraken.out" 
		kreport="${DIR}/classification/kraken2_fc_${DB}/${BASE}_filtered.kreport"
		export KRAKEN2_DB_PATH=/media/zbl/Storage/db/kraken2_standard
        
        # Check if mapping i/o files already exist
        if ! check_file "$kraken2_output" "kraken2 output $kraken2_output" ; then
            continue
        fi
        
        # Execute the command for paired-end mapping
        kraken2 --db "${DB}" --quick --threads 8 --report-zero-counts \
        --use-names --output "$kraken2_output" \
        --report "${BASE}.kreport" "$file"
        
        # Move Kraken2 Reports to folder
		echo "Copying Kraken2 Reports to folder" $(date -u)
		mv "/mnt/d/Guiuan/${BASE}.kreport" "$kreport"
        
        # Check for successful mapping i/o files
        if ! check_success "$kreport" "kraken2 report file $kreport"; then
            continue
        fi

        # Add verbose option to track progress
        echo "Processed alignment of $BASE with $DB." $(date -u)
done

# Run Bracken to estimate abundance
## Loops through all input files matching the pattern. 
for krakenreport in "${DIR}/classification/kraken2_fc_${DB}/"*.kreport; do
	# Extract BASE from the R1 filename
    BASE=$(basename "$krakenreport" | cut -d '_' -f 1,2)

		echo -e "\nAttempting to perform estimation of abundance.\n" $(date -u)

		kmer_distr= "/mnt/d/Guiuan/db/kraken2_standard/standard8/database150mers.kmer_distrib"
		kraken2_input="${DIR}/classification/kraken2_fc_${DB}/${BASE}_fc_${DB}.kreport"
		bracken_output= "${DIR}/classification/bracken_fc_${DB}/${BASE}_fc_${DB}.bracken"
		export KRAKEN2_DB_PATH="${DIR}/db/${DB}"
		
        # Check if mapping i/o files already exist
        if ! check_file "$bracken_output" "bracken output $bracken_output" ; then
            continue
        fi
                
		python "/mnt/d/Guiuan/Bracken-master/src/est_abundance.py" -i "$kraken2_input" \
		-k "$kmer_distr" -o "$bracken_output"
		
		# Check for successful mapping i/o files
        if ! check_success "$bracken_output" "bracken output $bracken_output"; then
            continue
        fi

done

echo -e "\nClassification pipeline finished." $(date -u)
echo "Elapsed time: $SECONDS seconds"
echo -e "Thank you for your patience.\n"
