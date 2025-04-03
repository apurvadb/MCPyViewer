# MCPyViewer 
This repository contains tools for visualizing contigs at viral integration sites using the results from searcHPV. The provided tools enable clear and informative representation of integration events, aiding in the interpretation of sequencing data.

If you utilize MCPyViewer in your analysis, please acknowledge our work by citing the following manuscript: 
Add citation

# searcHPV
A viral integration point detection tool for targeted capture sequencing data. Details on how to run this tool are described here : https://github.com/WenjinGudaisy/SearcHPV

## Getting started

1. Required resources
    - Unix like environment
2. Download and install the required resources
    - Download conda >=22.9.0: [https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
    - Download the "environment.yaml" file under this repository
    - Create conda environment for MCPyViewer:

          conda env create -f [your_path]/environment.yaml
    This command will automatically set up all the third-party tools and packages required for MCPyViewer. The name of the environment is "MCPV_analysis_toolkit" (edit later as needed).

    You can check the packages and tools in this environment by:

            conda list -n MCPV_analysis_toolkit

   You can update the environment by:

            conda env update -f [your_path]/environment.yaml

 3. Usage
        - Activate the conda environment

            conda activate MCPV_analysis_toolkit

 # Generate MCPV link plots 
Pipeline for generating a link plot of viral integration breakpoints within the human and MCPyV genomes, along with a plot illustrating the distribution of microhomology degree at MCPyV integration breakpoints (as shown in Figure 3A,B).

Usage:

    ./MCPV_link_plot.sh -w {workdir} -i samples.txt -o intermediate_files/

Note: samples.txt is a text file containing a list of samples, each listed in a single column with column name "Sample".
```
# Example contents of samples.txt
Sample
Sample1
Sample2
Sample3
Sample4
```
Outputs are stored in a directory called "MCV_link_plots" 
    



