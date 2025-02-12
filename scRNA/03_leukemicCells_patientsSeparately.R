library(Seurat)
library(tidyverse)
library(harmony)

library(clusterProfiler)
library(msigdbr)

library(ggpubr)
library(patchwork)
library(viridis)
library(aplot)

kdmm_palette <- c('RUNX1::RUNX1T1 knockdown' = '#f41626', 'mismatch control'= '#2538a5')

cbPalette <- c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "#999999")

vanGalenPalette <- c('GMP'='#4B0082', 'B'='#FFA500', 'T'="#009E73", 'Plasma'="#C1BC00", 'ProB'="#A5531C", 'Mono'="#008080", 'NA (MSC)'="#222222", 'ProMono'="#FF4856", 'cDC'="#7F7D70", 
                     'NK'="#7B68EE", 'lateEry'="#FF5B00", 'HSC'="#0000FF", 'CTL'="#00FF7F", 'earlyEry'='#DC143C', 'Prog'="#AF08AF", 'pDC' = "#5A9DFF")

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

# 1. Load in data and visualize what's already there ===========================

## 1.1. Load in data (Seurat object with leukemic cells) -----------------------

##### Generated in 01_QC_healthyLeukemicDistinction.R. 
##### A similar Seurat object (just supplemented with what this script adds to that object) is available at https://doi.org/10.5281/zenodo.14578307.
seu_leukemic <- readRDS(paste0(wd, '300_seu_leukemic.rds'))

## 1.2. UMAP by condition (Figure 4C) ------------------------------------------
pdf(paste0(wd, "330_leukemicOnly_umap_byCond.pdf"), height = 3.5, width = 8)
p <- map2(seu_leukemic, names(seu_leukemic),
          ~ DimPlot(.x, 
                    reduction = "umap",
                    cols = kdmm_palette,
                    group.by = "condition") +
            theme_void() +
            labs(title = '') +
            plot_annotation(title = .y))
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'bottom')
dev.off()

## 1.3. Violin plots of RUNX1T1 and  RUNX1::RUNX1T1 target genes (Figure 4H) ----

### 1.3.1. Preparations
genes_of_interest <- c('RUNX1T1', 'CD34', 'NFE2')

##### find adjusted p-values for differential expression of genes of interest
padj <- lapply(seu_leukemic, function(seuObj){
  Idents(seuObj) <- 'condition'
  de <- FindMarkers(seuObj,
                    features = genes_of_interest,
                    logfc.threshold = 0,
                    min.pct = 0,
                    ident.1 = 'RUNX1::RUNX1T1 knockdown',
                    assay = 'RNA')
  padj <- de$p_val_adj %>% signif(1)
  names(padj) <- rownames(de)
  padj <- c(padj[genes_of_interest[1]], padj[genes_of_interest[2]], padj[genes_of_interest[3]])
  return(padj)
})

seu_leukemic <- lapply(seu_leukemic, function(seuObj){
  Idents(seuObj) <- 'orig.ident'
  return(seuObj)
})

### 1.3.2. Plot
p <- map(seu_leukemic,
         ~ VlnPlot(.x %>% subset(downsample = 2000), 
                   features = genes_of_interest,
                   split.by = 'condition',
                   split.plot = TRUE,
                   cols = c(kdmm_palette[2], kdmm_palette[1]),
                   y.max = c(NA, NULL + 0.5),
                   assay = 'RNA',
                   combine = FALSE))

p <- map2(p, padj, function(plots_sublist, padj_vector){
  plots_sublist_new <- pmap(list(plots_sublist, padj_vector, names(padj_vector)), 
                            function(plt, padj_value, gene_name){
                              plt_new <- plt +
                                ggtitle(NULL) +
                                theme(legend.position = 'none',
                                      axis.title.x = element_blank(),
                                      axis.text.x = element_blank()) +
                                labs(y = paste0('Expression level ', gene_name)) +
                                geom_text(x = 1,
                                          y = 1,
                                          label = paste0('padj = ', as.character(padj_value)))
                              return(plt_new)
                            })
  return(plots_sublist_new)  
})

patient_panels <- map(p, 
                      ~ ggarrange(plotlist = .x, ncol = 1, nrow = 3, common.legend = TRUE, legend = 'bottom'))           


pdf(paste0(wd, "340_leukemicOnly_vln_RUNX1-ETO_targetGenes.pdf"), height = 7, width = 7)
ggarrange(plotlist = patient_panels, ncol = 3, nrow = 1, common.legend = FALSE)
dev.off()

# 2. Differential expression analysis ==========================================

## 2.1. Find all DE genes between the 2 conditions -----------------------------
seu_leukemic <- lapply(seu_leukemic, function(obj){
  Idents(obj) <- 'condition'
  return(obj)
})

markers <- map(seu_leukemic,
               ~ FindMarkers(.x,
                             ident.1 = 'RUNX1::RUNX1T1 knockdown',
                             assay = 'RNA',
                             logfc.threshold = 0,
                             min.pct = 0))

saveRDS(markers, paste0(wd, '328_leukemic_allMarkers_KDvsMM.rds'))
##### markers <- readRDS(paste0(wd, '328_leukemic_allMarkers_KDvsMM.rds'))

# 3. Compare the DE genes between scRNAseq and bulk RNAseq =====================

## 3.1. Prepare the names of significantly and substantially expressed genes from the bulk RNAseq ----

### 3.1.1. Load in the DESeq2 output of the bulk RNAseq of the PDX
##### Generated in this study in bulkRNA/bulkRNAseq.R

res_bulk <- read_delim('repos/runx1eto-kd-aml/bulkRNA/out/020_deseq2_out.csv') %>%
  rename (gene = `...1`)

### 3.1.2. Filter the bulk RNAseq results by significance and log2FC
res_bulk_filt <- res_bulk %>%
  dplyr::filter(padj < 0.05 & !is.na(padj) & abs(log2FoldChange) > 2)

bulk_up <- res_bulk_filt %>%
  dplyr::filter(log2FoldChange > 0) 

bulk_dn <- res_bulk_filt %>%
  dplyr::filter(log2FoldChange < 0) 

### 3.2.3. Prepare the input for GSEA
gsea_in_geneset <- data.frame(
  'gs_name' = paste0('bulkRNAseq patient C (PDX), ', 
                     c(rep('upregulated', nrow(bulk_up)), 
                       rep('downregulated', nrow(bulk_dn)))),
  'gene_symbol' = c(bulk_up$gene, bulk_dn$gene)
)

## 3.2. Prepare the ranks of DE genes from scRNAseq ----------------------------
ranks <- lapply(markers, function(df){
  df_filt <- df %>% filter(pct.1 + pct.2 > 0.01) ##### filter out genes whose expression is low
  ranks <- df_filt$avg_log2FC
  names(ranks) <- rownames(df_filt)
  ranks <- ranks[!is.na(ranks)] %>% 
    sort(decreasing = TRUE) ##### sorting is required for ClusterProfiler
  return(ranks)
})

## 3.3. Run GSEA ---------------------------------------------------------------
gsea_res <- map(ranks, ~ GSEA(
  geneList = .x,
  maxGSSize = 2000,
  pvalueCutoff = 1,
  eps = 0,
  pAdjustMethod = 'BH',
  seed = TRUE,
  TERM2GENE = gsea_in_geneset
))

csv_names <- paste0('343_gsea_res_againstBulk_', names(seu_leukemic), '.csv')
map2(gsea_res, csv_names, ~ write.csv(.x, paste0(wd, .y)))

saveRDS(gsea_res, paste0(wd, '343_gsea_res_againstBulk.rds'))
gsea_res <- readRDS(paste0(wd, '343_gsea_res_againstBulk.rds'))

## 3.4. Vizualize GSEA results (Suppl. Figure 4F) ------------------------------
pdf(paste0(wd, '344_gseaPlots_againstBulk.pdf'), width = 10, height = 10)
plotlist <- lapply(names(gsea_res), function(sample_name){
  sublist <- lapply(paste0('bulkRNAseq patient C (PDX), ', c('upregulated', 'downregulated')), function(geneset_name){
    p <- enrichplot::gseaplot(
      gsea_res[[sample_name]],
      geneSetID = geneset_name,
      title = paste0(sample_name, ', ', geneset_name),
      color.line = '#21908c') %>%
      as.patchwork()
    p <- p + plot_annotation(subtitle = paste0(
      'NES = ', as.character(signif(gsea_res[[sample_name]]@result$NES[gsea_res[[sample_name]]@result$ID == geneset_name], digits = 3)),
      ', padj = ', as.character(signif(gsea_res[[sample_name]]@result$p.adjust[gsea_res[[sample_name]]@result$ID == geneset_name], digits = 3))
    ))
    return(p)
  })
  return(sublist)
})
print(plotlist[[1]][[1]])
print(plotlist[[1]][[2]])
print(plotlist[[2]][[1]])
print(plotlist[[2]][[2]])
print(plotlist[[3]][[1]])
print(plotlist[[3]][[2]])
dev.off()

# 4. LSC score =================================================================

## 4.1. Add score --------------------------------------------------------------
##### Here, we are not using the original LSC6 formula but Seurat's AddModuleScore, in order to account for the sparcity of the single cell data.
LSC_signature <- c("DNMT3B", "CD34", "ADGRG1", "SOCS2", "SPINK2", "FAM30A")

seu_leukemic <- map(seu_leukemic,
                    ~ AddModuleScore(.x,
                                     features = list(LSC_signature),
                                     assay = 'RNA',
                                     name = 'LSC6'))

seu_leukemic <- lapply(seu_leukemic, function(seuObj){
  seuObj@meta.data <- dplyr::rename(seuObj@meta.data, LSC6 = LSC61)
  return(seuObj)
})

## 4.2. Plot pediatric LSC6 module score ---------------------------------------

### 4.2.1. Violin (Figure 7F)

##### plots to just get the p-values from
p1 <- map(seu_leukemic,
          ~ VlnPlot(.x, 
                    features = 'LSC6',
                    group.by = 'condition',
                    cols = c(kdmm_palette[2], kdmm_palette[1]),
                    y.max = c(NA, NULL + 0.5)) +
            stat_compare_means(comparisons = list(c('mismatch control', 'RUNX1::RUNX1T1 knockdown')), 
                               method = 'wilcox.test',
                               label = "p.format") +
            ggtitle(NULL) +
            theme(legend.position = 'none',
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank()))

panel1 <- ggarrange(plotlist = p1, ncol = 3, nrow = 1, common.legend = FALSE)

##### plots that I want but onto which I fail to put the significance bar
p2 <- map(seu_leukemic,
          ~ VlnPlot(.x %>% subset(downsample = 2000), 
                    features = 'LSC6',
                    split.by = 'condition',
                    split.plot = TRUE,
                    cols = c(kdmm_palette[2], kdmm_palette[1])) +
            ggtitle(NULL) +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank()) +
            ylim(-0.33, 0.68))

panel2 <- ggarrange(plotlist = p2, ncol = 3, nrow = 1, common.legend = TRUE, legend = 'bottom')

pdf(paste0(wd, '360_leukemicOnly_violin_byCond_pediatricLSC6.pdf'), width = 6.5, height = 5)
ggarrange(panel1, panel2, nrow = 2, ncol = 1)
dev.off()

##### code for boxplot, should you want to replace the violin plot with a box plot

# pdf(paste0(wd, '360_leukemicOnly_box_byCond_pediatricLSC6.pdf'), width = 6.5, height = 2.5)
# p <- map2(seu_leukemic, names(seu_leukemic),
#           ~ ggplot(.x@meta.data, 
#                     aes(x = condition, 
#                       y = LSC61, 
#                       fill = condition, 
#                       colour = condition)) +
#             geom_boxplot() +
#             theme_minimal() +
#             scale_fill_manual(values = kdmm_palette) +
#             scale_color_manual(values = c('black', 'black')) +
#             labs(y = 'Pediatric LSC6 score') +
#             theme(legend.position = 'none',
#                   panel.grid.major = element_blank(), 
#                   panel.grid.minor = element_blank(),
#                   axis.line = element_line(colour = "black"),
#                   axis.ticks = element_line(colour = "black"),
#                   axis.text.x = element_text(colour = "black"),
#                   axis.text.y = element_text(colour = "black")) +
#             stat_compare_means(comparisons = list(c('mismatch control', 'RUNX1/RUNX1T1 knockdown')), 
#                                method = 'wilcox.test',
#                                label = "p.format") +
#             plot_annotation(title = .y)) 
# ggarrange(plotlist = p, ncol = 3, common.legend = FALSE, legend = 'none')
# dev.off()

# 5. GSEA for the Zheng cord blood HSC gene set ================================

## 5.1. Get the gene set -------------------------------------------------------
pathways <- msigdbr(species = 'Homo sapiens', category = 'C8')
goi <- filter(pathways, 
              gs_name == 'ZHENG_CORD_BLOOD_C6_HSC_MULTIPOTENT_PROGENITOR'|gs_name == 'HAY_BONE_MARROW_CD34_POS_HSC') %>%
  select(gs_name, gene_symbol)

## 5.2. Run GSEA (ranks can be reused from step 3) -----------------------------
gsea_res <- map(ranks, ~ GSEA(
  geneList = .x,
  maxGSSize = 2000,
  pvalueCutoff = 1,
  eps = 0,
  pAdjustMethod = 'BH',
  seed = TRUE,
  TERM2GENE = goi
))

saveRDS(gsea_res, paste0(wd, '362_gsea_res_ZhengHSC_HayHSC.rds'))
gsea_res <- readRDS(paste0(wd, '362_gsea_res_ZhengHSC_HayHSC.rds'))

## 5.3. Vizualize GSEA results (Suppl. Figure 6A) ------------------------------
pdf(paste0(wd, '364_gseaPlots_OlafsRanks_Zheng_Hay_HSC.pdf'), width = 10, height = 10)
plotlist <- lapply(names(gsea_res), function(sample_name){
  sublist <- lapply(paste0(unique(goi$gs_name)), function(geneset_name){
    p <- enrichplot::gseaplot(
      gsea_res[[sample_name]],
      geneSetID = geneset_name,
      title = paste0(sample_name, ', ', geneset_name),
      color.line = '#21908c') %>%
      as.patchwork()
    p <- p + plot_annotation(subtitle = paste0(
      'NES = ', as.character(signif(gsea_res[[sample_name]]@result$NES[gsea_res[[sample_name]]@result$ID == geneset_name], digits = 3)),
      ', padj = ', as.character(signif(gsea_res[[sample_name]]@result$p.adjust[gsea_res[[sample_name]]@result$ID == geneset_name], digits = 3))
    ))
    return(p)
  })
  return(sublist)
})
print(plotlist[[1]][[1]])
print(plotlist[[1]][[2]])
print(plotlist[[2]][[1]])
print(plotlist[[2]][[2]])
print(plotlist[[3]][[1]])
print(plotlist[[3]][[2]])
dev.off()

pdf(paste0(wd, '364_gseaPlots_OlafsRanks_Hay_HSC.pdf'), width = 10, height = 10)
plotlist <- lapply(names(gsea_res), function(sample_name){
    p <- enrichplot::gseaplot(
      gsea_res[[sample_name]],
      geneSetID = 'HAY_BONE_MARROW_CD34_POS_HSC',
      title = paste0(sample_name, ', ', 'HAY_BONE_MARROW_CD34_POS_HSC'),
      color.line = '#21908c') %>%
      as.patchwork()
    p <- p + plot_annotation(subtitle = paste0(
      'NES = ', as.character(signif(gsea_res[[sample_name]]@result$NES[gsea_res[[sample_name]]@result$ID == 'HAY_BONE_MARROW_CD34_POS_HSC'], digits = 3)),
      ', padj = ', as.character(signif(gsea_res[[sample_name]]@result$p.adjust[gsea_res[[sample_name]]@result$ID == 'HAY_BONE_MARROW_CD34_POS_HSC'], digits = 3))
    ))
    return(p)
})
print(plotlist[[1]])
print(plotlist[[2]])
print(plotlist[[3]])
dev.off()

# 6. Cell type projections according to Zeng et al (2023) reference ============

##### Reference paper: https://www.biorxiv.org/content/10.1101/2023.12.26.573390v1
##### Reference Seurat object provided by Andy Zeng

## 6.1. Load in the reference and look at it -----------------------------------
ref <- readRDS('references/Zeng_BoneMarrowMap_Annotated_Dataset.rds')

Idents(ref) <- 'Tissue'

pdf(paste0(wd, '370_zengRef_UMAP_byCelltype.pdf'), width = 12)
DimPlot(ref %>% subset(downsample = 99999),
        reduction = 'umap',
        group.by = 'CellType_Broad',
        label = TRUE,
        cols = zengPalette) +
  theme_void() +
  labs(title = '')
dev.off()

## 6.2. Project ----------------------------------------------------------------
seu_leukemic <- lapply(seu_leukemic, function(seuObj) {
  seuObj <- NormalizeData(seuObj)
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)
  DefaultAssay(seuObj) <- 'RNA'
  return(seuObj)
})

anchors <- map(seu_leukemic, 
               ~ FindTransferAnchors(reference = ref, query = .x, dims = 1:30, npcs = 30))

anc1 <- FindTransferAnchors(reference = ref, query = seu_leukemic[[1]], dims = 1:30, npcs = 30)
saveRDS(anc1, paste0(wd, '376_anchrs_zeng_ptA.rds'))

anc2 <- FindTransferAnchors(reference = ref, query = seu_leukemic[[2]], dims = 1:30, npcs = 30)
saveRDS(anc2, paste0(wd, '376_anchrs_zeng_ptB.rds'))

anc3 <- FindTransferAnchors(reference = ref, query = seu_leukemic[[3]], dims = 1:30, npcs = 30)
saveRDS(anc3, paste0(wd, '376_anchrs_zeng_ptC.rds'))

anc1 <- readRDS(paste0(wd, '376_anchrs_zeng_ptA.rds'))
anc2 <- readRDS(paste0(wd, '376_anchrs_zeng_ptB.rds'))
anc3 <- readRDS(paste0(wd, '376_anchrs_zeng_ptC.rds'))
anchors <- list(anc1, anc2, anc3)
names(anchors) <- names(seu_leukemic)

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

##### Save predictions and prediction confidence scores because they will be needed later in the integrated object
saveRDS(preds, paste0(wd, '378_zengPreds.rds'))
preds <- readRDS(paste0(wd, '378_zengPreds.rds'))

seu_leukemic <- lapply(names(seu_leukemic), function(name){
  seuObj <- AddMetaData(seu_leukemic[[name]], metadata = preds[[1]][[name]], col.name = 'pred.Zeng.celltype') %>%
    AddMetaData(metadata = preds[[2]][[name]], col.name = 'pred.Zeng.score')
  return(seuObj)
})
names(seu_leukemic) <- names(preds[[1]])

## 6.3. Visualize --------------------------------------------------------------

### 6.3.1. UMAP by predictions 
pdf(paste0(wd, '386_leukemicOnly_umap_predZeng.pdf'), height = 3.5, width = 7)
p <- map2(seu_leukemic, names(seu_leukemic),
          ~ DimPlot(.x,
                    reduction = "umap",
                    cols = zengPalette,
                    group.by = "pred.Zeng.celltype") +
            theme(axis.line = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = 'bottom') +
            labs(title = NULL) +
            plot_annotation(title = .y))
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'bottom')
dev.off()

### 6.3.2. UMAP by prediction confidence 
pdf(paste0(wd, '394_leukemicOnly_umap_predZeng_score.pdf'), height = 3.5, width = 7)
p <- map2(seu_leukemic, names(seu_leukemic),
          ~ FeaturePlot(.x,
                        features = 'pred.Zeng.score') +
            scale_color_viridis() +
            theme(axis.line = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = 'bottom') +
            labs(title = NULL) +
            plot_annotation(title = .y))
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'bottom')
dev.off()

### 6.3.3. Box of prediction confidence by condition 
pdf(paste0(wd, '402_leukemicOnly_box_byCond_ZengScore.pdf'), width = 5, height = 2.5)
p <- map2(seu_leukemic, names(seu_leukemic),
          ~ ggplot(.x@meta.data, 
                   aes(x = condition, 
                       y = pred.Zeng.score, 
                       fill = condition, 
                       colour = condition)) +
            geom_boxplot() +
            theme_minimal() +
            scale_fill_manual(values = kdmm_palette) +
            scale_color_manual(values = c('black', 'black')) +
            labs(y = 'Zeng score') +
            theme(legend.position = 'none',
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black")) +
            stat_compare_means(comparisons = list(c('mismatch control', 'RUNX1::RUNX1T1 knockdown')), 
                               method = 'wilcox.test',
                               label = "p.format") +
            plot_annotation(title = .y))
ggarrange(plotlist = p, ncol = 3, common.legend = FALSE, legend = 'none')
dev.off() 

# 7. Cell type annotations according to van Galen & Hovestadt (2019) reference ====

##### reference paper: https://doi.org/10.1016/j.cell.2019.01.031
##### Reference Seurat object by Trincado et al, Biorxiv 2022, https://www.biorxiv.org/content/10.1101/2022.03.02.482638v1
##### Reference Seurat object downloaded from https://figshare.com/ndownloader/files/33951884

## 7.1. Load in the reference --------------------------------------------------
load('references/reference_vanGalen_VelascoHernandez.Rdata')
ref <- UpdateSeuratObject(merge.object)

## 7.2. Re-integrate the reference using Harmony and look at it ----------------
##### We do so to make it more similar to the Zeng et al. reference object.
ref <- RunHarmony(ref, 
                  group.by.vars = "orig.ident", 
                  dims.use = 1:30, 
                  max_iter = 50)
ref <- RunUMAP(ref, 
               reduction = "harmony", 
               dims = 1:30)
ref <- FindNeighbors(ref, reduction = "harmony", dims = 1:30)

pdf(paste0(wd, '410_refVanGalen_harmony_umap_byCellType.pdf'), width = 3.5, height = 3.5)
DimPlot(ref,
        group.by = 'predictionRF',
        cols = vanGalenPalette) +
  theme_void() +
  labs(title = '')
dev.off()

## 7.3. Project ----------------------------------------------------------------  
seu_leukemic <- lapply(seu_leukemic, function(seuObj) {
  seuObj <- NormalizeData(seuObj)
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)
  DefaultAssay(seuObj) <- 'RNA'
  return(seuObj)
})

anchors <- map(seu_leukemic, 
               ~ FindTransferAnchors(reference = ref, query = .x, dims = 1:30))

predictions <- map(anchors,
                   ~ TransferData(anchorset = .x, refdata = ref$predictionRF, dims = 1:30))

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
##### Save predictions and prediction confidence scores because they will be needed later in the integrated object
saveRDS(preds, paste0(wd, '418_vanGalenPreds_harmony.rds'))
##### preds <- readRDS(paste0(wd, '418_vanGalenPreds_harmony.rds'))

seu_leukemic <- lapply(names(seu_leukemic), function(name){
  seuObj <- AddMetaData(seu_leukemic[[name]], metadata = preds[[1]][[name]], col.name = 'pred.vanGalen.celltype') %>%
    AddMetaData(metadata = preds[[2]][[name]], col.name = 'pred.vanGalen.score')
  return(seuObj)
})
names(seu_leukemic) <- names(preds[[1]])

## 7.4. Visualize --------------------------------------------------------------

### 7.4.1. UMAP by predictions 
pdf(paste0(wd, '426_leukemicOnly_umap_predVanGalen.pdf'), height = 3.5, width = 8)
p <- map2(seu_leukemic, names(seu_leukemic),
          ~ DimPlot(.x,
                    reduction = "umap",
                    cols = vanGalenPalette,
                    group.by = "pred.vanGalen.celltype") +
            theme(axis.line = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = 'bottom') +
            labs(title = NULL) +
            plot_annotation(title = .y))
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'bottom')
dev.off()

### 7.3.2. UMAP by prediction confidence 
pdf(paste0(wd, '434_leukemicOnly_umap_predVanGalen_score.pdf'), height = 3.5, width = 8)
p <- map2(seu_leukemic, names(seu_leukemic),
          ~ FeaturePlot(.x,
                        features = 'pred.vanGalen.score') +
            scale_color_viridis() +
            theme(axis.line = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = 'bottom') +
            labs(title = NULL) +
            plot_annotation(title = .y))
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'bottom')
dev.off()

### 7.3.3. Box of prediction confidence by condition 
pdf(paste0(wd, '442_leukemicOnly_box_byCond_vanGalen_score.pdf'), width = 5, height = 2.5)
p <- map2(seu_leukemic, names(seu_leukemic),
          ~ ggplot(.x@meta.data, 
                   aes(x = condition, 
                       y = pred.vanGalen.score, 
                       fill = condition, 
                       colour = condition)) +
            geom_boxplot() +
            theme_minimal() +
            scale_fill_manual(values = kdmm_palette) +
            scale_color_manual(values = c('black', 'black')) +
            labs(y = 'van Galen score') +
            theme(legend.position = 'none',
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  axis.text.x = element_text(colour = "black"),
                  axis.text.y = element_text(colour = "black")) +
            stat_compare_means(comparisons = list(c('mismatch control', 'RUNX1::RUNX1T1 knockdown')), 
                               method = 'wilcox.test',
                               label = "p.format") +
            plot_annotation(title = .y))
ggarrange(plotlist = p, ncol = 3, common.legend = FALSE, legend = 'none')
dev.off() 

# 8. Make a list of Seurat objects for uploading to Zenodo =====================
seu_save <- lapply(seu_leukemic, function(obj){
  obj[['HTO']] <- NULL
  obj[['RNA']]$scale.data <- NULL
  obj[['SCT']]$counts <- NULL
  obj@meta.data <- obj@meta.data %>%
    select(condition, Phase, LSC6, pred.Zeng.celltype, pred.Zeng.score, pred.vanGalen.celltype, pred.vanGalen.score)
  return(obj)
})
saveRDS(seu_save, paste0(wd, 'seu_ptABCseparately_leukemicCells.rds'))

# 98. Session info ==============================================================
sink(paste0(wd, '999_sessionInfo_3.txt'))
sessionInfo()
sink()

# 99. Clear the environment =====================================================
rm(list = ls())