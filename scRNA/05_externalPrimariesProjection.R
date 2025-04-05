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

# 1. Load in the Zeng et al reference ==========================================
##### Reference paper: https://www.biorxiv.org/content/10.1101/2023.12.26.573390v1
##### Reference Seurat object provided by Andy Zeng
ref <- readRDS('references/Zeng_BoneMarrowMap_Annotated_Dataset.rds')
Idents(ref) <- 'Tissue'

# 2. RUNX1::RUNX1T1 AMLs from Lambo et al 2023 =================================

## 2.1. Load in the data -------------------------------------------------------
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

## 2.2. Subset to malignant cells only -----------------------------------------
seu <- lapply(seu, function(seuObj){
  seuObj_subset <- subset(seuObj, Malignant == 'Malignant')
  return(seuObj_subset)
})

## 2.2. Project onto Andy Zeng's reference -------------------------------------
seu <- lapply(seu, function(seuObj) {
  seuObj <- NormalizeData(seuObj)
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)
  DefaultAssay(seuObj) <- 'RNA'
  return(seuObj)
})

anchors <- map(seu, 
               ~ FindTransferAnchors(reference = ref, query = .x, dims = 1:30, npcs = 30))

anc16 <- FindTransferAnchors(reference = ref, 
                            query = seu[['AML16']], 
                            dims = 1:30, 
                            npcs = 30,
                            reference.reduction = "pca")
saveRDS(anc16, paste0(wd, '700_anchrs_zeng_lambo_aml16.rds'))

anc15 <- FindTransferAnchors(reference = ref, 
                            query = seu[['AML15']], 
                            dims = 1:30, 
                            npcs = 30,
                            reference.reduction = "pca")
saveRDS(anc15, paste0(wd, '700_anchrs_zeng_lambo_aml15.rds'))

anc14 <- FindTransferAnchors(reference = ref, 
                             query = seu[['AML14']], 
                             dims = 1:30, 
                             npcs = 30,
                             reference.reduction = "pca")
saveRDS(anc14, paste0(wd, '700_anchrs_zeng_lambo_aml14.rds'))

anc13 <- FindTransferAnchors(reference = ref, 
                             query = seu[['AML13']], 
                             dims = 1:30, 
                             npcs = 30,
                             reference.reduction = "pca")
saveRDS(anc13, paste0(wd, '700_anchrs_zeng_lambo_aml13.rds'))

anc12 <- FindTransferAnchors(reference = ref, 
                             query = seu[['AML12']], 
                             dims = 1:30, 
                             npcs = 30,
                             reference.reduction = "pca")
saveRDS(anc12, paste0(wd, '700_anchrs_zeng_lambo_aml12.rds'))

anc16 <- readRDS(paste0(wd, '700_anchrs_zeng_lambo_aml16.rds'))
anc15 <- readRDS(paste0(wd, '700_anchrs_zeng_lambo_aml15.rds'))
anc14 <- readRDS(paste0(wd, '700_anchrs_zeng_lambo_aml14.rds'))
anc13 <- readRDS(paste0(wd, '700_anchrs_zeng_lambo_aml13.rds'))
anc12 <- readRDS(paste0(wd, '700_anchrs_zeng_lambo_aml12.rds'))

anchors <- list(anc16, anc15, anc14, anc13, anc12)
names(anchors) <- names(seu)

seu <- map2(anchors, seu,
                ~ MapQuery(anchorset = .x, 
                           reference = ref, 
                           query = .y, 
                           refdata = list(Zeng_celltype = 'CellType_Broad'),
                           reference.reduction = "pca", 
                           reduction.model = "umap"))

preds <- lapply(seu, function(seuObj){
  df <- seuObj@meta.data
  df$cell_label <- rownames(df)
  df <- df %>% select(cell_label, predicted.Zeng_celltype, predicted.Zeng_celltype.score)
  return(df)
})

saveRDS(preds, paste0(wd, '701_zeng_lamboPreds.rds'))
preds <- readRDS(paste0(wd, '701_zeng_lamboPreds.rds'))

seu <- map2(seu, preds, function(seuObj, prediction_df){
  seuObj@meta.data$cell_label <- rownames(seuObj@meta.data)
  seuObj@meta.data <- seuObj@meta.data %>%
    left_join(prediction_df, by = 'cell_label')
  rownames(seuObj@meta.data) <- seuObj$cell_label
  return(seuObj)
})

pred_embeddings <- lapply(seu, function(seuObj){
  pca <- Embeddings(seuObj[["ref.pca"]])
  umap <- Embeddings(seuObj[["ref.umap"]])
  out <- list('pca' = pca, 'umap' = umap)
  return(out)
})

saveRDS(pred_embeddings, paste0(wd, '702_zeng_lamboPred_embeddings.rds'))
pred_embeddings <- readRDS(paste0(wd, '702_zeng_lamboPred_embeddings.rds'))

seu <- map2(seu, pred_embeddings, function(seuObj, embList){
  seuObj[['ref.pca']] <- CreateDimReducObject(as.matrix(embList[['pca']]))
  seuObj[['ref.umap']] <- CreateDimReducObject(as.matrix(embList[['umap']]))
  return(seuObj)
})

## 2.3. Visualize --------------------------------------------------------------
ref_data <- Embeddings(subset(ref, downsample = 50000)[["umap"]]) %>%
  as.data.frame() 

query_data <- lapply(names(seu), function(sample_name){
  df <- Embeddings(seu[[sample_name]][['ref.umap']]) %>%
    as.data.frame() %>%
    mutate(patient = sample_name) %>%
    rename(UMAP_1 = refUMAP_1,
           UMAP_2 = refUMAP_2) %>%
    tibble::rownames_to_column(var = 'cell_label') %>%
    left_join(seu[[sample_name]]@meta.data %>% select(cell_label, predicted.Zeng_celltype.score),
              by = 'cell_label')
  return(df)
})
names(query_data) <- names(seu)

plotlist <- lapply(names(seu), function(sample_name){
  plts <- query_data[[sample_name]] %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(data = ref_data, color = '#E3E3E3', size = 0.05, alpha = 0.5) +
      ggpointdensity::geom_pointdensity(size = 0.2) +
      jcolors::scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1.75) +
      geom_density_2d(alpha = 0.4, color = 'black', h = 1.5, linewidth = 0.3) +
      theme_void() + 
      labs(title = 'Sample from Lambo et al 2023 projected onto reference from Zeng et al 2023',
           subtitle = sample_name) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 18), 
                     legend.position = 'none')
    return(plts)
  })

pdf(paste0(wd, '705_lambo_umapAndy_density.pdf'), height = 4.5, width = 5)
plotlist
dev.off()

# 98. Session info ==============================================================
sink(paste0(wd, '999_sessionInfo_3.txt'))
sessionInfo()
sink()