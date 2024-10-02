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
bm <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = "hsapiens_gene_ensembl")

#import data
geneTrack_df <- read_csv("geneTrack_MCV.csv")
annotation_df <- read_csv("annotation_all_new.csv")
annotation_df = annotation_df[,-c(1)]
annotation_df = unique(annotation_df)
insertions_df <- read_csv("inserts_all_new.csv")
insertions_df = insertions_df[,-c(1)]
insertions_df = unique(insertions_df)

 
for (target_sample in unique(geneTrack_df$sample)) {
  for (target_gene in unique(geneTrack_df$gene)) {
    chr = filter(geneTrack_df,sample == target_sample, gene == target_gene)$chr[1]
    start = filter(geneTrack_df,sample == target_sample, gene == target_gene)$start[1]
    end = filter(geneTrack_df,sample == target_sample, gene == target_gene)$end[1]
    if (!is.na(chr)) {
      #build gene
      print(target_gene)
      print(target_sample)
      if (grepl("ENST",target_gene)){
        biomTrack <- BiomartGeneRegionTrack(genome = "hg38",
                                            name = "Genes",transcript  = target_gene,
                                            biomart = bm,
                                            col.line = NULL, col = NULL,transcriptAnnotation = "symbol")
      }
      else { biomTrack <- BiomartGeneRegionTrack(genome = "hg38",
                                                 name = "Genes",symbol = target_gene,
                                                 biomart = bm,
                                                 col.line = NULL, col = NULL,transcriptAnnotation = "symbol")
        
      }
     
      annotation_gene = filter(annotation_df,sampleID == target_sample,gene == target_gene)
      aaTrack <- AnnotationTrack(start = annotation_gene$start, width =annotation_gene$length, 
                                 chromosome = chr,
                                 genome = "hg38", name = annotation_gene$sampleName[1],shape = 'box',
                                 col = NULL,id = annotation_gene$strand,
                                 groupAnnotation = 'id',
                                 just.group = 'above',cex.group = 0.5,fill = 'black')
      axisTrack <- GenomeAxisTrack(add53=TRUE,add35=TRUE,labelPos = 'below')
      ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
      ht1 <- HighlightTrack(trackList = list(axisTrack),
                            start = unique(annotation_gene$ins), width = 1,
                            chromosome = chr, col = 'blue',fill = 'blue',inBackground = FALSE)
      
      #save plot
      
      pdf(file = paste('plots_new',paste(target_gene,target_sample,'pdf',sep = "."),sep = '/'),width = 10,height = 4.5)
      plotTracks(list(ideoTrack,ht1,biomTrack,aaTrack),from = start,to = end) 
      dev.off()
      #dev.copy2pdf(file = 'KLF13.pdf',width = 10,height = 3)
      
      #rstudioapi::savePlotAsImage("test.eps",format = c("eps"),width=1000,height=300)
      
      
      ####################
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
      insertions_gene$hpvSite = insertions_gene$hpvSite/5465*10
      insertions_gene$hpvEnd = insertions_gene$hpvEnd/5465*10
      insertions_gene$col2 = 'pink'
      insertions_gene[grep('E',insertions_gene$hpvGene),]$col2 = 'lightsalmon1'
      insertions_gene[grep('L',insertions_gene$hpvGene),]$col2 = 'turquoise3'
      
      insertions_gene[insertions_gene$type == 'hum',]$y = insertions_gene[insertions_gene$type == 'hum',]$y+0.1
      #xend1 = x + width
      #xend2 = x - width
      #xend = c(xend1,xend2)
      #y = rep(2,nrow(insertions_gene))
      #yend = rep(2,nrow(insertions_gene))
      #data = 'contig'
      #col = hpvGenome$col1[1:8]
      #insertions = data.frame(x = x,xend = xend,size = size,type = type,y = y,
      #yend = yend,data = data,hjust = hjust,col = insertions_gene$col,
      #dirc = insertions_gene$geneDirection)
      
      
      ##plot mcpv genome
      hpvGenome = read_csv('mcvDf.csv')
      hpvGenome = hpvGenome[,-1]
      colnames(hpvGenome) = c('type','direction','x','xend')
      hpvGenome$x = hpvGenome$x/5465*10
      hpvGenome$xend = hpvGenome$xend/5465*10
      hjust = (hpvGenome$xend+hpvGenome$x)/2
      vjust = runif(18,min = -5,max = -1)
      y = 1.2
      yend= 1.1
      size = 2
      data = 'mcpvGenome'
      hpvGenome = data.frame(hpvGenome,y,yend,size,data)
      hpvGenome$hjust = hjust
      col1 = c('darkgoldenrod1','darkgoldenrod1','mediumslateblue','mediumslateblue','mediumslateblue')
      hpvGenome$col1 = col1
      #hpvGenome$col2 = 'grey'
      #hpvGenome[grep('E',hpvGenome$type),]$col2 = 'darkgoldenrod1'
      #hpvGenome[grep('L',hpvGenome$type),]$col2 = 'mediumslateblue'
      #hpvGenome$type = recode(hpvGenome$type,`E1BS: E1 binding site` = 'E1BS', `E2BS: E2 binding site` = 'E2BS',
                              #`E6*`='E6',`E1^E4` = 'E1',`E8^E2` = 'E8')
      #hpvGenome = hpvGenome[-c(10,11),]
      #hpvGenome = hpvGenome[c(-8),]
      
      
      ######link segment hpv
      
      hpvshape = select(filter(insertions_gene,type == 'mcpv'),x,xend,hpvSite,hpvEnd,name,ins)
      x1 = melt(select(hpvshape,x,name,ins),id = c("name",'ins'),value.name = 'x')
      x2 = melt(select(hpvshape,xend,name,ins),id = c("name",'ins'),value.name = 'x')
      x3 = melt(select(hpvshape,hpvSite,name,ins),id = c("name",'ins'),value.name = 'x')
      x4 = melt(select(hpvshape,hpvEnd,name,ins),id = c("name",'ins'),value.name = 'x')
      hpvshape_x = rbind(x1,x2,x3,x4)
      hpvshape$y = 1.8
      hpvshape$yend = 1.3
      y1 = melt(select(hpvshape,y,name,ins),id = c("name",'ins'),value.name = 'y')
      y2 = melt(select(hpvshape,y,name,ins),id = c("name",'ins'),value.name = 'y')
      y3 = melt(select(hpvshape,yend,name,ins),id = c("name",'ins'),value.name = 'y')
      y4 = melt(select(hpvshape,yend,name,ins),id = c("name",'ins'),value.name = 'y')
      hpvshape_y = rbind(y1,y2,y3,y4)
      hpvshape = cbind(hpvshape_x,select(hpvshape_y,y)) %>% select(name,ins,x,y)
      
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
        geom_segment(hpvGenome,mapping = aes(x = x,xend=xend,y = y, yend = y),color = hpvGenome$col1,size = 4,
                  lineend = 'butt', linejoin = 'mitre', arrow = arrow(length = unit(insertions_gene$unit, "cm"),type = 'closed'))+
        geom_rect(hpvGenome,mapping = aes(xmin = 0,xmax=10,ymin = min(y)+0.05, ymax = max(yend)+0.05),color = 'grey28',
                  fill = NA,size = 0.02) +
        geom_text_repel(data = hpvGenome,mapping = aes(x = hjust, y = yend,label = type,fontface = 2),
                        color = 'black',size=3,max.overlaps = 30,min.segment.length = 0,nudge_y = -0.1,direction = 'x')+ ylim(0,2.5) + 
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank(),legend.position = "none") +
        geom_rect(filter(insertions_gene,type == 'hum'),mapping = aes(xmin = x,xmax=xend,ymin = y+0.1, ymax = y+0.2,fill = as.factor(ins)),alpha = 0.3,color='NA')+
        geom_rect(filter(insertions_gene,type == 'mcpv'),mapping = aes(xmin = x,xmax=xend,ymin = y-0.15, ymax = 1.8,fill = as.factor(ins)),alpha = 0.3,color = 'NA')+
        geom_polygon(hpvshape,mapping = aes(x = x,y = y,group = name,fill = as.factor(ins)),alpha = 0.3) + 
        geom_polygon(humshape,mapping = aes(x = x,y = y,group = name,fill = as.factor(ins)),alpha = 0.3)
      
      ggsave(paste('plots_new',paste(target_gene,target_sample,'1','pdf',sep = "."),sep = '/'),device ='pdf',plot = f,width = 10,height = 6)
      
      
      
      
    }
    
  }
}


 #+geom_segment(data = filter(insertions_gene,start == 0),mapping = aes(x = xend, y = y - 0.1, xend = xend, yend = 1.7),
#color = filter(insertions_gene,start == 0)$col2)+
#geom_segment(data = filter(insertions_gene,start == 0),mapping = aes(x = xend, 
#                                                                        y = y + 0.2, xend = xend, 
#                                                                        yend = y + 0.3),
#                color = 'tan', linetype = 'dashed') +
#   geom_segment(data = filter(insertions_gene,start == 0),mapping = aes(x = xend, y = y + 0.3, xend = (ins-geneStart)/(geneEnd-geneStart)*10, yend = 2.5),
#                color = 'tan', linetype = 'dashed') +
#   geom_segment(data = filter(insertions_gene,start == 0),mapping = aes(x = xend, y = 1.7, xend = hpvSite, yend = 1.2),
#                color = filter(insertions_gene,start == 0)$col2) 
# 



#dev.copy2pdf(,width = 10,height = 6)




#generate colors
#n <- 18
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#col1 = sample(col_vector, n)

#######################


