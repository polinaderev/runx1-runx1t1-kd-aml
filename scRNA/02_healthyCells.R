library(Seurat)
library(tidyverse)

library(EnhancedVolcano)
library(patchwork)
library(ggpubr)
library(viridis)

##### Another color-blind-friendly palette:
cbPalette2 <- colorBlindness::paletteMartin
names(cbPalette2) <- NULL
cbPalette2 <- c(cbPalette2[3:9], cbPalette2[11:15], "#999999")

##### Custom palette for our experimental condtions:
kdmm_palette <- c('RUNX1::RUNX1T1 knockdown' = '#f41626', 'mismatch control'= '#2538a5')

zengPalette <- c('B' = "#FFA500", 'Plasma Cell' = "#C1BC00", 'Plasma_Cell' = "#C1BC00", 'Stromal' = "#222222", 'Early Erythroid' = "#FF0000",
                 'HSC MPP' = "#0000FF", 'CD4 Memory T' = "#00FA9A", 'CD8 Memory T' = "#008000", 'cDC' = "#7F7D70",
                 'Cycling Progenitor' = "#A2EF65", 'Early GMP' = "#4B0082", 'Early Lymphoid' = "#800000", 
                 'EoBasoMast Precursor' = "#00CED1", 'Late Erythroid' = "#FF5B00", 'Late GMP' = "#FF69B4",
                 'LMPP' = "#2A53AD", 'Megakaryocyte Precursor' = "#00FF00", "MEP" = '#AF08AF', 'Monocyte' = "#008080",
                 'Naive T' = "#809999", 'NK' = "#7B68EE", 'pDC' = "#5A9DFF", 'Pre-B' = "#6B8E23", 'Pro-B' = "#A5531C",
                 'Pro-Monocyte' = "#FF4856")

##### Directory for outputs (replace with your directory for outputs):
wd <- 'repos/runx1eto-kd-aml/scRNA/out/'

##### reproducibility
set.seed(42)

# 1. Load in data ==============================================================

##### Seurat object with healthy cells
seu_healthy <- readRDS(paste0(wd, '290_seu_healthy.rds'))

# 2. Perform differential expression analysis ==================================

## 2.1 Find all markers between conditions (not included in the paper) ---------
seu_healthy <- lapply(seu_healthy, function(obj){
  Idents(obj) <- 'condition'
  return(obj)
})

markers <- map(seu_healthy,
               ~ FindMarkers(.x,
                             ident.1 = 'RUNX1::RUNX1T1 knockdown',
                             assay = 'RNA',
                             logfc.threshold = 0,
                             min.pct = 0))

saveRDS(markers, paste0(wd, '310_healthy_allMarkers_KDvsMM.rds'))
##### markers <- readRDS(paste0(wd, '310_healthy_allMarkers_KDvsMM.rds'))

## 2.2. Vizualize --------------------------------------------------------------

### 2.2.1. UMAP by condition (Figure 4C)
p <- map2(seu_healthy, names(seu_healthy),
          ~ DimPlot(.x, 
                    reduction = "umap",
                    cols = kdmm_palette,
                    group.by = "condition") +
            theme_void() +
            labs(title = .y))
pdf(paste0(wd, '315_healthy_umap_byCond.pdf'), height = 4)
ggarrange(plotlist = p, ncol = 2, common.legend = TRUE, legend = 'bottom')
dev.off()

### 2.2.2. Volcano plot (Suppl. Figure 4E)
markers <- map(markers, ~ dplyr::select(.x, avg_log2FC, p_val_adj))

pdf(paste0(wd, '320_healthy_volcano.pdf'))
map2(markers, names(markers),
     ~ EnhancedVolcano(.x,
                       lab = paste0("italic('", rownames(.x), "')"),
                       parseLabels = TRUE,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       drawConnectors = TRUE,
                       title = NULL,
                       subtitle = .y,
                       pointSize = c(ifelse(abs(.x$avg_log2FC) > 1 & .x$p_val_adj < 0.05, 3, 1)),
                       col = c('#444444','#444444', cbPalette2[7], cbPalette2[2]),
                       colAlpha = 1))
dev.off()

# 3. Perform projection onto the reference dataset =============================

##### Zeng et al (2025) reference 
##### reference paper: https://doi.org/10.1158/2643-3230.BCD-24-0342
##### reference Seurat object provided by Andy Zeng

## 3.1. Load in the reference --------------------------------------------------
ref <- readRDS('references/Zeng_BoneMarrowMap_Annotated_Dataset.rds')

## 3.2. Project ----------------------------------------------------------------
seu_healthy <- lapply(seu_healthy, function(seuObj) {
  seuObj <- NormalizeData(seuObj)
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)
  DefaultAssay(seuObj) <- 'RNA'
  return(seuObj)
})

anchors <- map(seu_healthy, 
               ~ FindTransferAnchors(reference = ref, query = .x, dims = 1:30, npcs = 30))

predictions <- map(anchors,
                   ~ TransferData(anchorset = .x, refdata = ref$CellType_Broad, dims = 1:30))

Idents(ref) <- 'CellType_Broad'

prediction.celltype <- lapply(names(predictions), function(name){
  pred <- factor(predictions[[name]]$predicted.id, levels = levels(ref))
  names(pred) <- rownames(predictions[[name]])
  return(pred)
})
names(prediction.celltype) <- names(predictions)

prediction.score <- lapply(names(predictions), function(name){
  pred <- predictions[[name]]$prediction.score.max
  names(pred) <- rownames(predictions[[name]])
  return(pred)
})
names(prediction.score) <- names(predictions)

preds <- list(prediction.celltype, prediction.score)

saveRDS(preds, paste0(wd, '322_healthy_zengPreds.rds'))
##### preds <- readRDS(paste0(wd, '322_healthy_zengPreds.rds'))

seu_healthy <- lapply(names(seu_healthy), function(name){
  seuObj <- AddMetaData(seu_healthy[[name]], metadata = preds[[1]][[name]], col.name = 'pred.Zeng.celltype') %>%
    AddMetaData(metadata = preds[[2]][[name]], col.name = 'pred.Zeng.score')
  return(seuObj)
})
names(seu_healthy) <- names(preds[[1]])

## 3.3. Visualize --------------------------------------------------------------

### 3.3.1. UMAP by prediction (Figure 4D)
pdf(paste0(wd, '323_healthy_ptA_umap_predZeng.pdf'), height = 6, width = 9)
DimPlot(seu_healthy[[1]],
        reduction = "umap",
        cols = zengPalette,
        group.by = "pred.Zeng.celltype",
        split.by = 'condition') &
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom') &
  labs(title = NULL)
dev.off()

### 3.3.2. UMAP by prediction confidence (Figure 4E)
pdf(paste0(wd, '324_healthy_umap_predZeng_score.pdf'), height = 5, width = 10)
map2(seu_healthy, names(seu_healthy),
     ~FeaturePlot(.x,
                  features = 'pred.Zeng.score',
                  split.by = 'condition') &
       scale_color_viridis(limits = c(min(.x@meta.data[['pred.Zeng.score']]), 
                                      max(.x@meta.data[['pred.Zeng.score']]))) &
       theme(axis.line = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks = element_blank(),
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             legend.position = 'bottom') &
       labs(subtitle = .y))
dev.off()

### 4.3.3. Box of prediction confidence by condition (Figure 4F, Suppl. Figure 4F)
pdf(paste0(wd, '325_healthy_box_byCond_ZengScore.pdf'), width = 2.5, height = 2.5)
map2(seu_healthy, names(seu_healthy),
     ~ ggplot(.x@meta.data, 
              aes(x = condition, 
                  y = pred.Zeng.score, 
                  fill = condition, 
                  colour = condition)) +
       geom_violin() +
       geom_boxplot(width = 0.2) +
       theme_minimal() +
       scale_fill_manual(values = kdmm_palette) +
       scale_color_manual(values = c('black', 'black')) +
       labs(y = 'Zeng score',
            subtitle = .y) +
       theme(legend.position = 'none',
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "black"),
             axis.ticks = element_line(colour = "black"),
             axis.text.x = element_text(colour = "black"),
             axis.text.y = element_text(colour = "black")) +
       stat_compare_means(comparisons = list(c('mismatch control', 'RUNX1::RUNX1T1 knockdown')), 
                          method = 'wilcox.test',
                          label = "p.format"))
dev.off()

# 5. Make a list of Seurat objects for uploading to Zenodo =====================
seu_save <- lapply(seu_healthy, function(obj){
  obj[['HTO']] <- NULL
  obj[['RNA']]$scale.data <- NULL
  obj[['SCT']]$counts <- NULL
  obj@meta.data <- obj@meta.data %>%
    select(condition, pred.Zeng.celltype, pred.Zeng.score)
  return(obj)
})
saveRDS(seu_save, paste0(wd, 'seu_ptABseparately_healthyCells.rds'))
##### Available at https://doi.org/10.5281/zenodo.14578307, folder "scRNA"

# 99. Session info ==============================================================
sink(paste0(wd, '999_sessionInfo_2.txt'))
sessionInfo()
sink()