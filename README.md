# MCPyViewer 
This repository contains tools for visualizing assembled contigs at viral-host integration sites using the results from searcHPV. SearcHPV was originally developed to assemble contigs around HPV-host integration sites [1] and later expanded to facilitate detection and contig assembly of multiple viruses that integrate into the human genome, including Merkel Cell Polyoma Virus (MCPyV) [2]. The provided tools enable clear and informative representation of integration events, aiding in the interpretation of sequencing data.

If you utilize MCPyViewer in your analysis, please acknowledge our work by citing the following manuscript: 

**Genomic Signatures of Poor Prognosis in Merkel Cell Carcinoma: A Single-Institution Prospective Study**

# Installation

    git clone https://github.com/apurvadb/MCPyViewer.git

# Prerequisite <h3>Running the "SearcHPV" pipeline</h3> 
searcHPV is a **viral integration detection** tool developed by our lab, designed to work with custom viral reference genomes. Prior to using MCPyViewer, please execute the searcHPV pipeline to generate viral integration breakpoint data, using the MCPyV genome (NC_010277.2) as the reference file. For step-by-step guidance on running the searcHPV pipeline, refer to the detailed instructions available here : 
```
https://github.com/WenjinGudaisy/SearcHPV
``` 
Once the pipeline has been succesfully executed, the resulting breakpoint data can be used for visualization. 

Please ensure that the resulting files are stored in the following folder structure in your working directory :

   
    {your_work_dir}/searcHPV_results/Sample_{your_sample_name}/

    

This is required for proper execution of the pipelines. 

## Getting started

1. Required resources
    - Unix like environment
2. Download and install the required resources
    - Download conda >=22.9.0: [https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
    - Download the "environment.yaml" file under this repository
    - Create conda environment for MCPyViewer:

          conda env create -f [your_path]/environment.yaml
    This command will automatically set up all the third-party tools and packages required for MCPyViewer. The name of the environment is "MCPyViewer".

    You can check the packages and tools in this environment by:

            conda list -n MCPyViewer

   You can update the environment by:

            conda env update -f [your_path]/environment.yaml

 3. Usage
        - Activate the conda environment

            conda activate MCPyViewer

All available pipelines for this toolkit are located in the "**scripts/**" folder. 

 # <h3>Generate MCPyV link plots</h3>
Pipeline for generating a link plot of viral integration breakpoints within the human and MCPyV genomes, along with a plot illustrating the distribution of degree of microhomology at MCPyV integration breakpoints (as shown in Figure 3A,B).

Usage:

    ./MCPV_link_plot.sh -w {workdir} -i samples.txt -o intermediate_files/

Note: The expected format for samples.txt is a text file with sample names listed in a single column, under the header "Sample".
```
# Example of samples.txt
Sample
Sample1
Sample2
Sample3
Sample4
```
Outputs are stored in a directory called "**MCPyV_link_plots**" 

Example output plots:
1. Link plot of viral integration breakpoints within the human and MCPyV genomes, with each line indicating the position of distinct integration breakpoints colored by MCPyV genes that the breakpoint fell within.

<p align="center">
    <img src="https://github.com/user-attachments/assets/519b6a0e-4539-4796-8a11-1ced7dad69d5" alt="Description" width="650">
</p> 

2. Distribution of degree of microhomology at breakpoints of MCPyV integrations. Number of overlapped base pairs of human and MCPyV segments at each breakpoint were calculated to represent microhomology. The
number of gapped base pairs at each breakpoint was calculated as negative score, with clean breaks denoted as zero-base pair overlapped.

<p align="center">
    <img src="https://github.com/user-attachments/assets/af7e6043-9689-448e-a434-9077c56d8b75" alt="Description" width="500">
</p> 


# <h3>Generate MCPyV gene model plots</h3>
Pipeline for generating MCPyV integration gene model (as shown in Fig 3C and Supplementary Figure S4).

Usage :

    ./MCPV_geneModel.sh -w <workdir> -t <path_to_transcript_gtf> -e <path_to_exon_gtf> -r <path_to_reference_fa> -d <path_to_output_dir> -f <ideogram_hg38_file> <sample1> <sample2> ...

Note: 
The "ideogram_hg38_file (ideogram_hg38_data.txt)", "transcript_gtf (Homo_sapiens.GRCh38.105.transcript.gtf)", "exon_gtf (Homo_sapiens.GRCh38.105.exon.gtf)", "reference_fa (hg_mcv.fa)" files are available in the "**data/**" folder.

If you have not already indexed your merged Human+MCPV reference, please do so by following these commands:

```
#activate MCPyViewer conda environment first to make sure you are using the correct versions of tools
ref = '[path_of_your_reference_file]'
bwa index {ref}
samtools faidx {ref}
java -jar $PICARDLIB/picard.jar CreateSequenceDictionary R={ref} O={ref.replace('.fa','.dict')
```
We are utilizing the MCPyV reference from NCBI, which can be accessed at https://www.ncbi.nlm.nih.gov/nuccore/NC_010277.2/. An example reference file available for use can be found at data/hg_mcv.fa.

Please note that this entire section can be executed as a batch script on a cluster, as it utilizes BWA for sequence alignment, which may demand additional memory and resources.

Example :

```
#!/bin/bash
#SBATCH --job-name=geneModel
#SBATCH --mail-user=xxx@xxx.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40gb
#SBATCH --time=100:00:00
#SBATCH --account=XXXXX
#SBATCH --partition=standard
#SBATCH --output=geneModel.log
#SBATCH --error=geneModel.err

conda activate MCPyViewer
./MCPV_geneModel.sh -w <workdir> -t <path_to_transcript_gtf> -e <path_to_exon_gtf> -r <path_to_reference_fa> -d <path_to_output_dir> -f <ideogram_hg38_file> <sample1> <sample2> ...    

```

If you encounter an "Ensembl website unresponsive" message, please rerun the pipeline. This issue may arise due to temporary connectivity problems with the Ensembl website.

Outputs are stored in a directory called "**MCPyV_geneModel_plots**"

Example output :

Representative MCPyV integration events in a tumor.

<p align="center">
    <img src="Images/geneModel.png" alt="Description" width="800">
</p> 

## <h4>Test run on publicly available data [3].</h4> 
Available here: https://www.ncbi.nlm.nih.gov/sra/ERX4366251

Running from the "scripts/" folder

Usage :

```
1. ./MCPV_link_plot.sh -w test/ -i test/test_sample.txt -o test/<your_output_dir>

2. ./MCPV_geneModel.sh -w test/ -t <path to>/data/Homo_sapiens.GRCh38.105.transcript.gtf -e <path to>/data/Homo_sapiens.GRCh38.105.exon.gtf -r <path to>/data/hg_mcv.fa -d test/{your_output_dir} -f <path to>/data/ideogram_hg38_data.txt ERR4425693

```

test_sample.txt, MCPV integration analysis results are located in the test folder under "searcHPV_results/Sample_ERR4425693". 
Outputs will be available in "MCPyV_link_plots" and "MCPyV_geneModel_plots" folders. 

# References 

[1]: Pinatti, Lisa M et al. “SearcHPV: A novel approach to identify and assemble human papillomavirus-host genomic integration events in cancer.” Cancer vol. 127,19 (2021): 3531-3540. doi:10.1002/cncr.33691
     https://pubmed.ncbi.nlm.nih.gov/34160069/        

[2]: Genomic Signatures of Poor Prognosis in Merkel Cell Carcinoma: A Single-Institution Prospective Study

[3]: Czech-Sioli, Manja et al. “High-resolution analysis of Merkel Cell Polyomavirus in Merkel Cell Carcinoma reveals distinct integration patterns and suggests NHEJ and MMBIR as underlying mechanisms.” PLoS pathogens vol. 16,8 e1008562. 24 Aug. 2020, 
     doi:10.1371/journal.ppat.1008562
