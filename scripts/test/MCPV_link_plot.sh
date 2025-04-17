#!/bin/bash

# Usage function to display help
usage() {
    echo "Usage: $0 -i <integration_data> -s <samples_file> -o <output_dir>"
    exit 1
}

# Parse command-line arguments
while getopts ":i:s:o:" opt; do
    case ${opt} in
        i )
            integration_data=$OPTARG
            ;;
        s )
            samples_file=$OPTARG
            ;;
        o )
            output_dir=$OPTARG
            ;;
        \? )
            echo "Invalid option: -$OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Check if all arguments are provided
if [ -z "$integration_data" ] || [ -z "$samples_file" ] || [ -z "$output_dir" ]; then
    echo "Missing required argument(s)"
    usage
fi

# Run the Python script
python input_to_link_plots.py -i "$integration_data" -s "$samples_file" -o "$output_dir"
if [ $? -ne 0 ]; then
    echo "Python script failed"
    exit 1
fi

# Run the R script
Rscript MCPV_plots.R -d "$output_dir/" -o "$output_dir/MCPyV_link_plots/"
if [ $? -ne 0 ]; then
    echo "R script failed"
    exit 1
fi

echo "Both scripts executed successfully"
