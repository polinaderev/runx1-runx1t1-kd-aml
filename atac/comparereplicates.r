#Normalise as CPM in excel - x*1000000/sum(tagcount_allpeaks)
#Save tag counts for all peaks or distal peaks (>1.5kb away from TSS)

R

setwd("Z:/LauraSwart/Peaks")
library(pheatmap)
all <- read.csv("Normalisedpeaks_all.csv",header=TRUE)
distal <- read.csv("Normalisedpeaks_distal.csv",header=TRUE)
all <- log2(all+0.1)
distal <-  log2(distal+0.1)
corall <- cor(all,method="pearson")
cordistal <- cor(distal,method="pearson")
pheatmap(corall,border_color=F,show_colnames = F)
pheatmap(cordistal,border_color=F,show_colnames = F)
