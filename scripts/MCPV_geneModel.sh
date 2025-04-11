#!/bin/bash

# Helper function to display usage information
usage() {
    echo "Usage: $0 -w workdir -t transcript_gtf -e exon_gtf -r ref_fa -d output_dir -f file_to_process sample1 sample2 ... sampleN"
    exit 1
}

# Initialize variMCV_geneModelles
WORKDIR=""
TRANSCRIPT_GTF=""
EXON_GTF=""
REF_FA=""
OUTPUT_DIR=""
FILE_TO_PROCESS=""
SAMPLES=()

# Print received arguments for debugging
echo "Received arguments: $@"

# Parse command-line arguments
while getopts ":w:t:e:r:d:f:" opt; do
  case $opt in
    w) WORKDIR="$OPTARG";;
    t) TRANSCRIPT_GTF="$OPTARG";;
    e) EXON_GTF="$OPTARG";;
    r) REF_FA="$OPTARG";;
    d) OUTPUT_DIR="$OPTARG";;
    f) FILE_TO_PROCESS="$OPTARG";;
    \?) echo "Invalid option: -$OPTARG" >&2; usage;;
    :) echo "Option -$OPTARG requires an argument." >&2; usage;;
  esac
done

# Shift to move past the processed options
shift $((OPTIND - 1))

# Remaining arguments should be samples
SAMPLES=("$@")

# Print variMCV_geneModelle values for debugging
echo "WORKDIR: $WORKDIR"
echo "TRANSCRIPT_GTF: $TRANSCRIPT_GTF"
echo "EXON_GTF: $EXON_GTF"
echo "REF_FA: $REF_FA"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "FILE_TO_PROCESS: $FILE_TO_PROCESS"
echo "SAMPLES: ${SAMPLES[@]}"

# Check if all required arguments are provided
if [ -z "$WORKDIR" ] || [ -z "$TRANSCRIPT_GTF" ] || [ -z "$EXON_GTF" ] || [ -z "$REF_FA" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$FILE_TO_PROCESS" ] || [ "${#SAMPLES[@]}" -eq 0 ]; then
    echo "Error: Missing required arguments."
    usage
else
    echo "All required arguments are provided."
fi

# Run the Python script
python input_to_geneModel.py -w "${WORKDIR}" -g "${TRANSCRIPT_GTF}" -e "${EXON_GTF}" -r "${REF_FA}" -o "${OUTPUT_DIR}" -s "${SAMPLES[@]}"

# Check if the first command was successful
if [ $? -ne 0 ]; then
    echo "Error: input_to_geneModel.py failed."
    exit 1
fi

# Run the R script
Rscript MCPV_geneModel.R -d "${OUTPUT_DIR}" -o MCPyV_geneModel_plots/ -f "${FILE_TO_PROCESS}"

# Check if the second command was successful
if [ $? -ne 0 ]; then
    echo "Error: MCPV_geneModel.R failed."
    exit 1
fi

echo "Scripts executed successfully."
