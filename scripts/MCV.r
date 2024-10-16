library(circlize)
library(stringr)
library(RColorBrewer)
library(dplyr)
library(readr)
library(ggplot2)
library(reshape2)


#primary link plot

pdf('link_plot.pdf',width = 13,height = 11)
circos.clear()
orig_track <- read.delim("MCV_track.bed",header = FALSE,stringsAsFactors = FALSE)
colnames(orig_track) = c('gene','start','end')
hpvTrack = orig_track[1:8,]
humTrack = filter(orig_track,gene != 'chrX'&gene != 'chrY')
humTrack  = humTrack[9:nrow(humTrack),]
humTrack$newGene =gsub("chr","",humTrack$gene)
humTrack = humTrack[order(as.numeric(humTrack$newGene)),]
temp = rbind(hpvTrack[order(hpvTrack$start),],humTrack[,1:3])
track = rbind(temp,filter(orig_track,gene == 'chrX'|gene == 'chrY'))


hpvSite = read.delim("MCVIntegration.bed",header = FALSE,stringsAsFactors = FALSE)
geneSite = read.delim("genomeIntegration_MCV.bed",header = FALSE,stringsAsFactors = FALSE)
colnames(hpvSite) = c('gene','start','end')
colnames(geneSite) = c('gene','start','end')
label = read.delim("label.bed",header = FALSE,stringsAsFactors = FALSE)
colnames(label) = c('chr','start','end','value1')


col_12_sample = brewer.pal(5, 'Set3')
col2 = rep("lightBlue",24)
col = c(col_12_sample[1:3],col_12_sample[1],col_12_sample[4],col_12_sample[1],col_12_sample[5],col_12_sample[4],col2)

circos.genomicInitialize(track,labels.cex = 0.5,track.height = 0.2,plotType = NULL)
circos.track(ylim = c(0, 3), 
             bg.col = col, 
             bg.border = NA, track.height = 0.1, panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               if(str_detect(chr,"chr") == FALSE){
                 label = paste(chr,sep = " ")
               }
               else{
                 label = chr
               }
               
               circos.text(mean(xlim),mean(ylim)+2,label , cex = 1, col = "Black",
                           facing = "clockwise", niceFacing = TRUE,adj = c(-0.1))
             })


hpvGenes = unique(track$gene[1:7])

linkColor = rep(0,nrow(hpvSite))
i = 1
for(each in hpvSite$gene){
  j = 1
  for(eachGene in hpvGenes){
    if(each == eachGene){
      linkColor[i] = c(col_12_sample[1:3],col_12_sample[1],col_12_sample[4],col_12_sample[1],col_12_sample[5],col_12_sample[4])[j]
    }
    j = j + 1
  }
  i = i + 1
}

circos.genomicLink(hpvSite,geneSite, col = linkColor, border = NA,lwd = 2)
dev.off()
#primary integration fell into gene
twelveSamplePrimaryGene <- read_csv("/Volumes/GoogleDrive/My Drive/research/hpv fusion/project/12_samples/twelveSamplePrimaryGene.csv")
colnames(twelveSamplePrimaryGene) = c('gene','count')
ggplot(twelveSamplePrimaryGene) + geom_bar(aes(gene,count),stat = 'identity',fill = 'orange',color = 'white') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10),text = element_text(size=20)) 
                                                                                                        

#hpv count
twelveSamplePrimaryHPV <- read_csv("MCV_gene.csv")
colnames(twelveSamplePrimaryHPV) = c('gene','observed','expected')
newTwelveSamplePrimaryHPV = melt(twelveSamplePrimaryHPV)
colnames(newTwelveSamplePrimaryHPV) = c('gene','legend','count')
virus_gene <- ggplot(newTwelveSamplePrimaryHPV) + geom_bar(aes(gene,count,fill = legend),stat = 'identity',position = 'dodge') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20), legend.title = element_blank(), axis.text.x = element_text(angle = 45,vjust = 0.5)) 
ggsave('MCV_gene.png')
#chi square test
twelveSamplePrimaryHPV$prob = twelveSamplePrimaryHPV$expected/sum(twelveSamplePrimaryHPV$expected)

#overall
ob = chisq.test(twelveSamplePrimaryHPV$observed,p = twelveSamplePrimaryHPV$prob)

#each
gene = rep(0,length(twelveSamplePrimaryHPV))
pvalue = rep(0,length(twelveSamplePrimaryHPV))
for (i in 1:nrow(twelveSamplePrimaryHPV)){
  pvalue[i] = chisq.test(c(twelveSamplePrimaryHPV$observed[i],sum(twelveSamplePrimaryHPV$observed)-twelveSamplePrimaryHPV$observed[i]),
             p = c(twelveSamplePrimaryHPV$prob[i],1-twelveSamplePrimaryHPV$prob[i]))$p.value
  gene[i] = twelveSamplePrimaryHPV$gene[i]
}

gene
pvalue

######
###microhomology
scores = read_csv("score.csv")
ggplot(NULL) + geom_histogram(aes(scores),stat = 'count',fill = 'orange') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20), legend.title = element_blank(), axis.text.x = element_text(vjust = 0.5,face = 'bold'),axis.text.y = element_text(face = 'bold'),
        axis.title.x = element_text(face = 'bold'),axis.title.y = element_text(face = 'bold')) + scale_x_continuous(breaks=c(seq(-40,20,10))) +
  xlab('#bp of overlap') + ylab('frequency')
ggsave('microhonology.pdf',width = 10,height = 10,dpi = 300,units = c('in'),device = 'pdf')


