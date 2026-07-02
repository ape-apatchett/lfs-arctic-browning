#!/bin/bash

# Set the input and output directories
#Latnja browning study metagenomics
input_dir="/xxx/raw"
output_dir="/xxx/Metagenomics/"


# Set the path to Trimmomatic
trimmomatic_path="/usr/share/java/trimmomatic-0.39.jar"

#make unpaired directory
mkdir -p "unpaired_dir"

# Iterate over the input files
#for input_file in "$input_dir"/*R1.fastq; do
  # Extract the file name without extension
#  file_name=$(basename "$input_file" _R1.fastq)

for input_file in "$input_dir"/*_R1.fastq.gz; do
  # Extract the file name without extension
  file_name=$(basename "$input_file" _R1.fastq.gz)

  # Set the input file paths
 # r1_input="$input_dir/${file_name}_R1.fastq"
 # r2_input="$input_dir/${file_name}_R2.fastq"

  # Set the corresponding input file paths
  r1_input="$input_file"
  r2_input="$input_dir/${file_name}_R2.fastq.gz"

  # Set the output file paths
  r1_output="$output_dir/${file_name}_R1_trimmed.fastq"
  r1_unpaired="$output_dir/${file_name}_R1_unpaired.fastq"
  r2_output="$output_dir/${file_name}_R2_trimmed.fastq"
  r2_unpaired="$output_dir/${file_name}_R2_unpaired.fastq"

  # Run Trimmomatic to trim and filter the reads including adapter file
  #java -jar "$trimmomatic_path" PE -threads 4 -phred33 "$r1_input" "$r2_input" "$r1_output" "$r1_unpaired" "$r2_output" "$r2_unpaired" \
#ILLUMINACLIP:/xxx/trimmomatic/adapters/allAdapter.fas:2:30:10 \
 #LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

#In the ILLUMINACLIP command above 2 = seed mismatches, allows up to 2 mismatches when matching the adapter sequence
#30 = Palindrome clip threshold, threshold score for paired-end reads
#10 =  simple clip threshold, threshold score for single-end reads

# Run Trimmomatic to trim and filter the reads (stricter quality and length settings, will likely remove more low-quality reads and short reads)
  java -jar "$trimmomatic_path" PE -threads 4 -phred33 "$r1_input" "$r2_input" "$r1_output" "$r1_unpaired" "$r2_output" "$r2_unpaired" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50

# Move untrimmed files to the "unpaired" directory
  mv "$r1_unpaired" "$unpaired_dir/"
  mv "$r2_unpaired" "$unpaired_dir/"

done


