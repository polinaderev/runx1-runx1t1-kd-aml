#Get differential peaks in R

library(dplyr)
library(DESeq2)
counts <- read.delim("Counts.tsv",comment.char = "#", header=TRUE)
colnames(counts) <- colnames(counts) %>% gsub(pattern="alignments.",replacement="")
colnames(counts) <- colnames(counts) %>% gsub(pattern="_rmdup.bam",replacement="")
countsfordeseq <- counts[,7:ncol(counts)]
coldata <- read.delim("coldata.txt",sep="\t",header=TRUE,row.names=1)
condFactor<- factor(coldata$condFactor, levels=(c('siRE','siMM')))
dds<- DESeqDataSetFromMatrix(countData = countsfordeseq, design = ~ condFactor, colData=coldata)
dds<- DESeq(dds)
results<- as.data.frame(results(dds, contrast = c('condFactor', 'siRE', 'siMM')))
coords<- counts[,2:4]
final<-cbind(coords,results)
write.table(final,"Differential_peaks.tsv",col.names=TRUE,sep="\t")

#Get ordered by FC, significantly gained, lost peaks in excel, save text files
