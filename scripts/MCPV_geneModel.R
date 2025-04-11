library(Gviz)
library(biomaRt)
library(lattice)
library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(stringr)
library(reshape2)
library(tidyr)
library(cowplot)
library(optparse)

#rm(list = ls())

# Setup command-line options
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default=getwd(),
              help="Directory containing all required input files", metavar="character"),
  make_option(c("-f", "--ideogram_file"), type="character", default=NULL,
              help="Path to ideogram hg38 file", metavar="character"),
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
geneTrack_file <- file.path(opt$directory, "geneTrack_MCV.csv")
annotation_file <- file.path(opt$directory, "annotation.csv")
inserts_file <- file.path(opt$directory, "inserts.csv")
ideogram_file <- file.path(opt$ideogram_file)
mcvDf_file <- file.path(opt$directory, "mcvDf.csv")


#httr::set_config(httr::config(ssl_verifypeer = 0, ssl_verifyhost = 2))
#setwd("/gpfs/accounts/chadbren_root/chadbren/apurvadb/MCPV_searcHPV/MCPyV Scripts")
bm <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
geneTrack_df <- read_csv(geneTrack_file)
annotation_df <- unique(read_csv(annotation_file))
insertions_df <- unique(read_csv(inserts_file))


ideo_data = read.table(ideogram_file, sep = "\t", header = FALSE)
colnames(ideo_data) = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
#ideo_data$chrom = gsub("^chr", "", ideo_data$chrom)

for (target_sample in unique(geneTrack_df$sample)) {
  for (target_gene in unique(geneTrack_df$gene)) {
    gene_data <- filter(geneTrack_df, sample == target_sample, gene == target_gene)
    chr <- gene_data$chr[1]
    start <- gene_data$start[1]
    end <- gene_data$end[1]

    if (!is.na(chr)) {
      print(target_gene)
      print(target_sample)

      biomTrack <- BiomartGeneRegionTrack(genome = "hg38", 
                                          name = "Genes",
                                          biomart = bm, 
                                          transcriptAnnotation = "symbol", 
                                          chromosome = gene_data$chr, start = gene_data$start, end = gene_data$end)

      annotation_gene <- filter(annotation_df, sampleID == target_sample, gene == target_gene)
      aaTrack <- AnnotationTrack(start = annotation_gene$start, width = annotation_gene$length,
                                 chromosome = annotation_gene$chr, genome = "hg38",
                                 name = annotation_gene$sampleName[1], 
                                 shape = 'box',
                                 col = NULL, 
                                 id = annotation_gene$strand,
                                 groupAnnotation = 'id',
                                 just.group = 'above', cex.group = 0.5, fill = 'black'
                                 )
      
      axisTrack <- GenomeAxisTrack(add53 = TRUE, add35 = TRUE, labelPos = 'below')
      ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = annotation_gene$chr, bands = ideo_data)
      ht1 <- HighlightTrack(trackList = list(axisTrack),
                            start = unique(annotation_gene$ins), width = 1,
                            chromosome = annotation_gene$chr, col = 'blue', fill = 'blue', inBackground = FALSE)

      pdf(file = paste0(opt$output, target_gene, '.', target_sample, '.pdf'), width = 11, height = 4.5)
      plotTracks(list(ideoTrack, axisTrack, biomTrack, aaTrack, ht1), from = start, to = end, col.line = NULL, col = NULL)
      dev.off()
      plot1_grob <- grid.grabExpr(
        plotTracks(list(ideoTrack, axisTrack, biomTrack, aaTrack, ht1), from = start, to = end, col.line = NULL, col = NULL)
      ) 


      insertions_gene = filter(insertions_df,sample == target_sample,gene == target_gene)
      insertions_gene$geneStart = min(insertions_gene$geneStart)
      insertions_gene$geneEnd = max(insertions_gene$geneEnd)
      insertions_gene[is.na(insertions_gene$feature),]$feature = 'Intergenic'
      insertions_gene = unique(insertions_gene)
      insertions_gene = insertions_gene[with(insertions_gene,order(sample,name,start)),]
      insertions_gene$col = 'red'
      insertions_gene[insertions_gene$type == 'hum',]$col = 'black'
      insertions_gene$size = 3
      #insertions_gene[insertions_gene$type == 'mcpv',]$size = 2
      insertions_gene$unit= 0.1
      insertions_gene[insertions_gene$feature == 'intron',]$unit = 0.1
      start = insertions_gene$start
      width =10/length(unique(insertions_gene$name))/3
      x = seq(0,10,10/length(unique(insertions_gene$name)))
      x = x[1:length(unique(insertions_gene$name))]
      x = x[match(insertions_gene$name, unique(insertions_gene$name))]
      insertions_gene$x = x + insertions_gene$start/insertions_gene$fullLength*width
      insertions_gene$xend = x + insertions_gene$end/insertions_gene$fullLength*width
      vjust = rep(c(2,-2),nrow(insertions_gene)/2 + 1)
      insertions_gene$name
      insertions_gene$vjust = vjust[1:nrow(insertions_gene)]
      insertions_gene$hjust = 0
      insertions_gene$y = 2
      insertions_gene$mcpvSite = insertions_gene$mcpvSite/5465*10
      insertions_gene$mcpvEnd = insertions_gene$mcpvEnd/5465*10
      insertions_gene$col2 = 'pink'
      insertions_gene[grep('E',insertions_gene$mcpvGene),]$col2 = 'lightsalmon1'
      insertions_gene[grep('L',insertions_gene$mcpvGene),]$col2 = 'turquoise3'
      
      insertions_gene[insertions_gene$type == 'hum',]$y = insertions_gene[insertions_gene$type == 'hum',]$y+0.1
      
      
      ##plot mcpv genome
      mcpvGenome = read_csv(mcvDf_file)
      #mcpvGenome = mcpvGenome[,-1]
      colnames(mcpvGenome) = c('type','direction','x','xend')
      mcpvGenome$x = mcpvGenome$x/5465*10
      mcpvGenome$xend = mcpvGenome$xend/5465*10
      hjust = (mcpvGenome$xend+mcpvGenome$x)/2
      vjust = runif(18,min = -5,max = -1)
      y = 1.2
      yend= 1.1
      size = 2
      data = 'mcpvGenome'
      mcpvGenome = data.frame(mcpvGenome,y,yend,size,data)
      mcpvGenome$hjust = hjust
      col1 = c('darkgoldenrod1','darkgoldenrod1','mediumslateblue','mediumslateblue','mediumslateblue')
      mcpvGenome$col1 = col1
      
      
      
      ######link segment mcpv
      
      mcpvshape = select(filter(insertions_gene,type == 'mcpv'),x,xend,mcpvSite,mcpvEnd,name,ins)
      x1 = melt(select(mcpvshape,x,name,ins),id = c("name",'ins'),value.name = 'x')
      x2 = melt(select(mcpvshape,xend,name,ins),id = c("name",'ins'),value.name = 'x')
      x3 = melt(select(mcpvshape,mcpvSite,name,ins),id = c("name",'ins'),value.name = 'x')
      x4 = melt(select(mcpvshape,mcpvEnd,name,ins),id = c("name",'ins'),value.name = 'x')
      mcpvshape_x = rbind(x1,x2,x3,x4)
      mcpvshape$y = 1.8
      mcpvshape$yend = 1.3
      y1 = melt(select(mcpvshape,y,name,ins),id = c("name",'ins'),value.name = 'y')
      y2 = melt(select(mcpvshape,y,name,ins),id = c("name",'ins'),value.name = 'y')
      y3 = melt(select(mcpvshape,yend,name,ins),id = c("name",'ins'),value.name = 'y')
      y4 = melt(select(mcpvshape,yend,name,ins),id = c("name",'ins'),value.name = 'y')
      mcpvshape_y = rbind(y1,y2,y3,y4)
      mcpvshape = cbind(mcpvshape_x,select(mcpvshape_y,y)) %>% select(name,ins,x,y)
      
      ######link segment hum
      humshape = select(filter(insertions_gene,type == 'hum'),x,xend,name,ins,geneStart,geneEnd,length,y)
      humshape$x3 = (humshape$ins-humshape$geneStart)/(humshape$geneEnd-humshape$geneStart)*10
      humshape$x4 = humshape$x3 + humshape$length/(humshape$geneEnd-humshape$geneStart)*10
      x1 = melt(select(humshape,x,name,ins),id = c("name",'ins'),value.name = 'x')
      x2 = melt(select(humshape,xend,name,ins),id = c("name",'ins'),value.name = 'x')
      x3 = melt(select(humshape,x3,name,ins),id = c("name",'ins'),value.name = 'x')
      x4 = melt(select(humshape,x4,name,ins),id = c("name",'ins'),value.name = 'x')
      humshape_x = rbind(x1,x2,x4,x3)
      humshape$y = humshape$y + 0.2
      humshape$yend = 2.5
      y1 = melt(select(humshape,y,name,ins),id = c("name",'ins'),value.name = 'y')
      y2 = melt(select(humshape,y,name,ins),id = c("name",'ins'),value.name = 'y')
      y3 = melt(select(humshape,yend,name,ins),id = c("name",'ins'),value.name = 'y')
      y4 = melt(select(humshape,yend,name,ins),id = c("name",'ins'),value.name = 'y')
      humshape_y = rbind(y1,y2,y3,y4)
      humshape = cbind(humshape_x,select(humshape_y,y)) %>% select(name,ins,x,y)
      
      
      ##plot contigs
       
       
      f = ggplot() +
        geom_segment(data = insertions_gene,mapping = aes(x = x, y = y, xend = xend, yend = y),color = insertions_gene$col,
                     lineend = 'butt', linejoin = 'mitre',
                     size = insertions_gene$size,arrow = arrow(length = unit(insertions_gene$unit, "cm"),type = 'closed'))+
        geom_text(data = insertions_gene,mapping = aes(x = x, y = y,label = feature,fontface = 2),
                  hjust = 0,vjust = insertions_gene$vjust,color = 'black',size=3) +
        geom_segment(mcpvGenome,mapping = aes(x = x,xend=xend,y = y, yend = y),color = mcpvGenome$col1,size = 4,
                  lineend = 'butt', linejoin = 'mitre', arrow = arrow(length = unit(insertions_gene$unit, "cm"),type = 'closed'))+
        geom_rect(mcpvGenome,mapping = aes(xmin = 0,xmax=10,ymin = min(y)+0.05, ymax = max(yend)+0.05),color = 'grey28',
                  fill = NA,size = 0.02) +
        geom_text_repel(data = mcpvGenome,mapping = aes(x = hjust, y = yend,label = type,fontface = 2),
                        color = 'black',size=3,max.overlaps = 30,min.segment.length = 0,nudge_y = -0.1,direction = 'x')+ ylim(0,2.5) + 
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank(),legend.position = "none") +
        geom_rect(filter(insertions_gene,type == 'hum'),mapping = aes(xmin = x,xmax=xend,ymin = y+0.1, ymax = y+0.2,fill = as.factor(ins)),alpha = 0.3,color='NA')+
        geom_rect(filter(insertions_gene,type == 'mcpv'),mapping = aes(xmin = x,xmax=xend,ymin = y-0.15, ymax = 1.8,fill = as.factor(ins)),alpha = 0.3,color = 'NA')+
        geom_polygon(mcpvshape,mapping = aes(x = x,y = y,group = name,fill = as.factor(ins)),alpha = 0.3) + 
        geom_polygon(humshape,mapping = aes(x = x,y = y,group = name,fill = as.factor(ins)),alpha = 0.3) +
        theme(plot.margin = margin(0, 0, 0, 20))
      
      ggsave(paste(opt$output,paste(target_gene,target_sample,'1','pdf',sep = "."),sep = '/'),device ='pdf',plot = f,width = 10,height = 6)
      combined_plot = plot_grid(plot1_grob, f, ncol = 1, rel_heights = c(0.55, 0.55))
      ggsave(paste0(opt$output, target_gene, '.', target_sample, '.combined.pdf'),
             combined_plot, width = 11, height = 11)
      
    }
    
  }
}


