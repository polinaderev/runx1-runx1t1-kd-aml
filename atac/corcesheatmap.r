#Make heatmap of Corces cell types vs RL048

library(pheatmap)
library(dplyr)
corces <- read.delim("Distalsites_Corcesunion_annotated.bed", sep="\t",header=TRUE)
corces <- corces[,20:ncol(corces)]
logc <- log2(corces+0.1)
logc[is.na(logc)] <- 0
cor<-cor(logc,method="pearson")
cor <- cor[7:14,]
cor <- cor[,1:6]
rownames(cor) <- rownames(cor) %>% gsub(pattern="X.rds.projects.b.boniferc.exphaem.Public_Data.Corces_et_al_ATACseq.CellType_peaks.",replacement="")
rownames(cor) <- rownames(cor) %>% gsub(pattern="_treat_pileup.bdg.bedGraph.avg.over.200.bp",replacement="")
colnames(cor) <- colnames(cor) %>% gsub(pattern="_treat_pileup.bdg.bedGraph.avg.over.200.bp",replacement="")
colnames(cor) <- colnames(cor) %>% gsub(pattern="Peaks.",replacement="")
pheatmap(cor,border_color=F)
