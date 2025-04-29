#!/bin/bash

SECONDS=0

echo -e "#####################################################################################

  ZBL Metagenomics Pre-processing Pipeline

  This bash script is optimized for the preprocessing of raw metagenomic 
	reads from Illumina sequencing (PGC).

  Usage:
	Enter PATH to directory where raw reads are stored. Indicate one folder up. 
	
  Example:
	PE reads: /media/zbl/Storage/mgs(/raw)

  Output:
  	1. Cleaned, trimmed, and error corrected PE reads in ~/path/clean
	2. Bowtie 2 output files and filtered reads in ~/path/filter_mouse


  ©Z.B. Lara, 2024

######################################################################################"

# Function to check file existence and proceed accordingly
check_file() {
    if [ -f "$1" ]; then
        echo "Skipping $2: File $1 already exists."
        return 1
    else
        echo "Error: File $1 does not exist."
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

#Assign variables
read -p "Enter PATH to files: " DIR
	if [ -d "$DIR" ]; then
	echo "The directory exists."
	else echo "The directory does not exist."
	fi &&
	
	echo "Creating output directories."
	mkdir "${DIR}/clean"
	mkdir -p "${DIR}/filter_human"

#Activate conda environment
echo -e "\nActivating metagenomics conda environment."
source ~/miniconda3/bin/activate mgs
conda env list &&

#QC
echo -e "Attempting to process files.\n"

# Loop through all R1 files matching the pattern *_R1_001.fastq.gz
for R1_file in "${DIR}/raw/"*_1.fq.gz; do
    
    # Extract BASE from the R1 filename
    BASE=$(basename "$R1_file" _1.fq.gz)

    # Construct the corresponding R2 filename
    R2_file="${DIR}/raw/${BASE}_2.fq.gz"

    # Check if the corresponding R2 file exists
    if [ -e "$R2_file" ]; then
    	
    	# QUALITY CONTROL
    	echo -e "\nAttempting to perform qc.\n"
       
        # Define the output filenames for QC
        qc_output_R1="${DIR}/clean/${BASE}_R1.fastq.gz"
        qc_output_R2="${DIR}/clean/${BASE}_R2.fastq.gz"
        json="${DIR}/clean/${BASE}_fastp_report.json"
        html_report="${DIR}/clean/${BASE}_fastp_report.html"
        
        # Check if QC output files already exist
		if [[ -f "$qc_output_R1" && -f "$qc_output_R2" ]]; then
			echo -e "\nQC outputs already exist for $BASE. Skipping fastp QC step.\n" $(date -u)
		else
			# Execute the command for paired-end qc
			fastp -i "$R1_file" -o "$qc_output_R1" -I "$R2_file" -O "$qc_output_R2" \
			--correction -p -w 8 --dedup -R "$json" -h "$html_report"
		fi

		# MAPPING
		echo -e "\nAttempting to perform filtering of human genome from reads.\n" $(date -u)

		# Define the output filenames for mapping      
		sam_file="${DIR}/filter_human/${BASE}_filter_mouse_bt2.sam"
		log_file="${DIR}/filter_human/${BASE}_filter.stats.log"       
		filtered_1="${DIR}/filter_human/${BASE}_filtered_R1.fastq.gz"
		filtered_2="${DIR}/filter_human/${BASE}_filtered_R2.fastq.gz"
  		export BOWTIE2_INDEXES="${DIR}/db"

		# Execute the command for paired-end mapping
		(bowtie2 -x GRCh38 -1 "$qc_output_R1" -2 "$qc_output_R2" \
		--fast-local -p 8 -t --no-unal --un-conc-gz "${BASE}_filtered" \
		-S "$sam_file") 2> "$log_file" &&

		echo -e "\nCopying mapped reads to folder.\n" $(date -u)

		mv "${DIR}/${BASE}_filtered.1" "$filtered_1" & 
		mv "${DIR}/${BASE}_filtered.2" "$filtered_2" &&

		# Check if mapping i/o files already exist
		if ! check_success "$filtered_2" "filtered output file $filtered_2"; then
			continue
		fi        	
				
		echo -e "\n$BASE pre-processing finished" $(date -u)

# add verbose option to track progress
echo "Processed paired-end files for $BASE." $(date -u)
        
        #add verbose option to track progress
        echo "Processed paired-end files for $BASE."  $(date -u)
    else
        echo "Corresponding R2 file not found for $R1_file."
    fi
done

echo -e "\nPipeline finished without errors" $(date -u)
echo "Elapsed time: $SECONDS seconds"
echo -e "Thank you for your patience.\n"
