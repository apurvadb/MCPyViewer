# MCPyViewer 
This repo includes scripts to visualize the contig at virus integration from SearchMCPV - edit this part 

If you use the MCPyViewer visualization tool, please cite our manuscript:
Add citation

# searchMCPV
A MCPV integration point detection tool for targeted capture sequencing data. Details on how to run this tool are described here : https://github.com/WenjinGudaisy/SearcHPV![image](https://github.com/user-attachments/assets/07543b07-c423-41e4-81a2-6452ab564939)
- update to refelct MCPV genome (Ask Wenjin)

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
 Pipeline to generate link plot of viral integration breakpoints within the human and MCPyV genomes and plot showing the distribution of degree of microhomology at breakpoints of MCPyV integrations (shown in Figure 3 (A,B))

    ./MCPV_link_plot.sh -w {workdir} -i samples.txt -o intermediate_files/

Example samples.txt file :
    Sample
    tumor1
    tumor2
    tumor3
    .
    .
    tumorN

     
         Outputs are stored in a directory called "MCV_link_plots" 
    



