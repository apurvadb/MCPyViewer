#!/bin/bash

# Help message
usage() {
    echo "Usage: $0 -w <workdir> -i <input_file> -o <output_dir>"
    exit 1
}

# Parse command-line arguments
while getopts "w:i:o:" opt; do
    case "$opt" in
        w) workdir="$OPTARG" ;;
        i) input_file="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$workdir" || -z "$input_file" || -z "$output_dir" ]]; then
    usage
fi

# Run Python script
python input_to_link_plots.py -w "$workdir" -i "$input_file" -o "$output_dir"

# Run R script
Rscript MCPV_plots.R -d "$output_dir/" -o "MCPyV_link_plots/"


