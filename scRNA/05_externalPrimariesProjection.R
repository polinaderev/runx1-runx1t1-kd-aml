library(Seurat)
library(Matrix)
library(tidyverse)

library(patchwork)
library(ggpubr)
library(viridis)

cbPalette <- c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "#999999")

zengPalette <- c('B' = "#FFA500", 'Plasma Cell' = "#C1BC00", 'Plasma_Cell' = "#C1BC00", 'Stromal' = "#222222", 'Early Erythroid' = "#FF0000",
                 'HSC MPP' = "#0000FF", 'CD4 Memory T' = "#00FA9A", 'CD8 Memory T' = "#008000", 'cDC' = "#7F7D70",
                 'Cycling Progenitor' = "#A2EF65", 'Early GMP' = "#4B0082", 'Early Lymphoid' = "#800000", 
                 'EoBasoMast Precursor' = "#00CED1", 'Late Erythroid' = "#FF5B00", 'Late GMP' = "#FF69B4",
                 'LMPP' = "#2A53AD", 'Megakaryocyte Precursor' = "#00FF00", "MEP" = '#AF08AF', 'Monocyte' = "#008080",
                 'Naive T' = "#809999", 'NK' = "#7B68EE", 'pDC' = "#5A9DFF", 'Pre-B' = "#6B8E23", 'Pro-B' = "#A5531C",
                 'Pro-Monocyte' = "#FF4856")

##### Directory for outputs (replace with your directory):
wd <- 'repos/runx1-runx1t1-kd-aml/scRNA/out/'

##### reproducibility
set.seed(42)

# 1. RUNX1::RUNX1T1 AMLs from Lambo et al 2023 =================================

## 1.1. Load in the data -------------------------------------------------------
##### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235063
ref_path <- 'references/Lambo2023/'

sample_names <- paste0(
  c('GSM7494257_AML16_',
    'GSM7494266_AML15_',
    'GSM7494314_AML14_',
    'GSM7494329_AML13_',
    'GSM7494326_AML12_'),
  'DX_processed_'
)

seu <- lapply(sample_names, function(sample_name){
  mtx <- readMM(paste0(ref_path, sample_name, 'matrix.mtx'))
  genes <- read.delim(paste0(ref_path, sample_name, "genes.tsv"),
                      header = FALSE, 
                      stringsAsFactors = FALSE)
  barcodes <- read.delim(paste0(ref_path, sample_name, "barcodes.tsv"),
                         header = FALSE, 
                         stringsAsFactors = FALSE)
  metadata <- read.delim(paste0(ref_path, sample_name, "metadata.tsv"), 
                         header = TRUE, 
                         row.names = 1)
  rownames(mtx) <- genes$V1  
  colnames(mtx) <- barcodes$V1
  seuObj <- CreateSeuratObject(counts = mtx, meta.data = metadata)
  return(seuObj)
})

names(seu) <- c('AML16', 'AML15', 'AML14', 'AML13', 'AML12')



# 98. Session info ==============================================================
sink(paste0(wd, '999_sessionInfo_3.txt'))
sessionInfo()
sink()