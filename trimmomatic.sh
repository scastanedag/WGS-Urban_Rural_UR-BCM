#!/bin/bash

# Sergio Castaneda
# Check if an input directory is provided
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input_directory>"
  exit 1
fi

# Variables
INPUT_DIR="$1"  # Input directory containing .fastq.gz files
THREADS=16  # Number of threads for Trimmomatic
ADAPTERS="TruSeq3-PE.fa"  # Path to adapter file, ensure it's available or specify the full path

# Check if the input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "The input directory $INPUT_DIR does not exist."
  exit 1
fi

# Create an output directory inside the input directory
OUTPUT_DIR="${INPUT_DIR}/trimmed_output"
mkdir -p "$OUTPUT_DIR"

# Iterate over all files ending with .1.fastq.gz in the input directory
for file1 in "$INPUT_DIR"/*.1.fastq.gz; do
  # Extract the base name of the file (without the suffix)
  base=$(basename "$file1" .1.fastq.gz)

  # Define the corresponding paired file
  file2="${INPUT_DIR}/${base}.2.fastq.gz"

  # Check if the paired file exists
  if [[ -f "$file2" ]]; then
    # Define the output file paths
    paired1="${OUTPUT_DIR}/${base}.1_paired_trim.fq.gz"
    unpaired1="${OUTPUT_DIR}/${base}.1_unpaired_trim.fq.gz"
    paired2="${OUTPUT_DIR}/${base}.2_paired_trim.fq.gz"
    unpaired2="${OUTPUT_DIR}/${base}.2_unpaired_trim.fq.gz"

    # Run Trimmomatic
    echo "Processing $file1 and $file2..."
    trimmomatic PE -threads $THREADS \
      "$file1" "$file2" \
      "$paired1" "$unpaired1" \
      "$paired2" "$unpaired2" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10:2:keepBothReads \
      MINLEN:150 AVGQUAL:20 TRAILING:20

    echo "Files processed and saved to $OUTPUT_DIR"
  else
    echo "The paired file $file2 does not exist. Skipping $file1."
  fi
done

echo "All files have been processed. Output saved in $OUTPUT_DIR."
