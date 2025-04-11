#!/usr/bin/env Rscript

# Function to ensure all packages are installed
#install_if_missing <- function(packages) {
#  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
#  if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
#}

# Specify the needed packages
required_packages <- c("circlize", "stringr", "RColorBrewer", "dplyr", "readr", "ggplot2", "reshape2", "optparse")

# Install missing packages
#install_if_missing(required_packages)

# Load the packages
lapply(required_packages, library, character.only = TRUE)

#library(circlize)
#library(stringr)
#library(RColorBrewer)
#library(dplyr)
#library(readr)
#library(ggplot2)
#library(reshape2)
#library(optparse)

# Setup command-line options
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default=getwd(),
              help="Directory containing all required input files", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="output",
              help="Directory to save output files [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Create output directory if it does not exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Construct full paths to input files
mcv_track_file <- file.path(opt$directory, "MCV_track.bed")
mcv_integration_file <- file.path(opt$directory, "MCVIntegration.bed")
genome_integration_file <- file.path(opt$directory, "genomeIntegration_MCV.bed")
label_file <- file.path(opt$directory, "label.bed")
mcv_gene_file <- file.path(opt$directory, "MCV_gene.csv")
score_file <- file.path(opt$directory, "score.csv")

# Output files
pdf_output <- file.path(opt$output, "link_plot.pdf")
mcv_gene_png <- file.path(opt$output, "MCV_gene.png")
microhomology_pdf <- file.path(opt$output, "microhomology.pdf")

# Read input files
orig_track <- read.delim(mcv_track_file, header = FALSE, stringsAsFactors = FALSE)
colnames(orig_track) <- c('gene','start','end')
mcvTrack <- orig_track[1:8,]
humTrack <- filter(orig_track, gene != 'chrX' & gene != 'chrY')
humTrack  <- humTrack[9:nrow(humTrack),]
humTrack$newGene <- gsub("chr", "", humTrack$gene)
humTrack <- humTrack[order(as.numeric(humTrack$newGene)),]
temp <- rbind(mcvTrack[order(mcvTrack$start),], humTrack[,1:3])
track <- rbind(temp, filter(orig_track, gene == 'chrX' | gene == 'chrY'))

mcvSite <- read.delim(mcv_integration_file, header = FALSE, stringsAsFactors = FALSE)
geneSite <- read.delim(genome_integration_file, header = FALSE, stringsAsFactors = FALSE)
colnames(mcvSite) <- c('gene', 'start', 'end')
#replace_dict = {'intergenic_1' : 'intergenic', 'intergenic_2' : 'intergenic', 'intergenic_3' : 'intergenic', 'large_T_1' : 'large_T', 'large_T_2' : 'large_T'}
#mcvSite['gene'] = mcvSite['gene'].replace(replace_dict)

colnames(geneSite) <- c('gene', 'start', 'end')
label <- read.delim(label_file, header = FALSE, stringsAsFactors = FALSE)
colnames(label) <- c('chr', 'start', 'end', 'value1')

col_12_sample <- brewer.pal(5, 'Set3')
col2 <- rep("lightBlue", 24)
col <- c(col_12_sample[1:3], col_12_sample[1], col_12_sample[4], col_12_sample[1],
          col_12_sample[5], col_12_sample[4], col2)

pdf(pdf_output, width = 13, height = 11)
circos.clear()

circos.genomicInitialize(track, labels.cex = 0.5, track.height = 0.2, plotType = NULL)
circos.track(ylim = c(0, 3), 
             bg.col = col, 
             bg.border = NA, track.height = 0.1, panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               if (str_detect(chr, "^chr")) {
                 label = chr  # Keep chromosomes as they are (chr1, chr2, etc.)
               } else {
                 # Extract the label from the 4th column
                 label = label[label[[1]] == chr, 4]
                 
                 # If no label found, default to chr name
                 if (length(label) == 0 || is.na(label)) {
                   label = chr
                 }
                 
                 # Ensure formatting as per requirement
                 label = paste(label, sep = " ")
               }      
	       #if (str_detect(chr, "chr") == FALSE) {
               #  label = paste(chr, sep = " ")
               #} else {
               #  label = chr
               #}
               circos.text(mean(xlim), mean(ylim) + 2, label, cex = 1, col = "Black",
                           facing = "clockwise", niceFacing = TRUE, adj = c(-0.1))
             })

mcvGenes <- unique(track$gene[1:7])
linkColor <- rep(0, nrow(mcvSite))
i <- 1
for(each in mcvSite$gene) {
  j <- 1
  for(eachGene in mcvGenes) {
    if(each == eachGene) {
      linkColor[i] = c(col_12_sample[1:3], col_12_sample[1], col_12_sample[4], col_12_sample[1],
                        col_12_sample[5], col_12_sample[4])[j]
    }
    j <- j + 1
  }
  i <- i + 1
}
circos.genomicLink(mcvSite, geneSite, col = linkColor, border = NA, lwd = 2)
dev.off()

# Reading MCV gene data
twelveSamplePrimaryMCV <- read_csv(mcv_gene_file)
colnames(twelveSamplePrimaryMCV) <- c('gene', 'observed', 'expected')
newTwelveSamplePrimaryMCV <- melt(twelveSamplePrimaryMCV)
colnames(newTwelveSamplePrimaryMCV) <- c('gene', 'legend', 'count')

# Generate and save MCV gene plot
virus_gene <- ggplot(newTwelveSamplePrimaryMCV) + 
  geom_bar(aes(gene, count, fill = legend), stat = 'identity', position = 'dodge') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size = 20), legend.title = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.5))
ggsave(mcv_gene_png)

# Chi square test
twelveSamplePrimaryMCV$prob <- twelveSamplePrimaryMCV$expected / sum(twelveSamplePrimaryMCV$expected)
#ob <- chisq.test(twelveSamplePrimaryMCV$observed, p = twelveSamplePrimaryMCV$prob)
twelveSamplePrimaryMCV$prob

# Chi square test for each gene
#gene <- rep(0, length(twelveSamplePrimaryMCV$gene))
#pvalue <- rep(0, length(twelveSamplePrimaryMCV$gene))
#for (i in 1:nrow(twelveSamplePrimaryMCV)) {
#  pvalue[i] <- chisq.test(
#    c(twelveSamplePrimaryMCV$observed[i], sum(twelveSamplePrimaryMCV$observed) - twelveSamplePrimaryMCV$observed[i]),
#    p = c(twelveSamplePrimaryMCV$prob[i], 1 - twelveSamplePrimaryMCV$prob[i])
#  )$p.value
#  gene[i] <- twelveSamplePrimaryMCV$gene[i]
#}
#print(gene)
#print(pvalue)

# Microhomology histogram
score_data <- read_csv(score_file)
ggplot(score_data, aes(x = Scores)) +
  geom_histogram(binwidth = 1, fill = "orange", color = "black") +  # Binwidth auto-adjusts for score range
  scale_x_continuous(breaks = seq(min(score_data$Scores), max(score_data$Scores), by = 2)) + # Dynamic axis scaling
  labs(
    #title = "Distribution of Scores",
    x = "#bp of overlap",
    y = "frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold") 
    #panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()
  )


ggsave(microhomology_pdf, width = 10, height = 10, dpi = 300, units = 'in', device = 'pdf')
