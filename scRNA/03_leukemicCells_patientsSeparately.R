library(Seurat)
library(tidyverse)

library(clusterProfiler)
library(msigdbr)

library(ggpubr)
library(patchwork)
library(viridis)
library(aplot)

kdmm_palette <- c('RUNX1::RUNX1T1 knockdown' = '#f41626', 'mismatch control'= '#2538a5')

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

# 1. Load in data and visualize what's already there ===========================

## 1.1. Load in data (Seurat object with leukemic cells) -----------------------

##### Generated in 01_QC_healthyLeukemicDistinction.R. 
##### A similar Seurat object (just supplemented with what this script adds to that object) is available at https://doi.org/10.5281/zenodo.14578307.
seu_leukemic <- readRDS(paste0(wd, '300_seu_leukemic.rds'))

## 1.2. UMAP by condition (Figure 4B) ------------------------------------------
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

## 1.3. Violin plots of RUNX1T1 and  RUNX1::RUNX1T1 target genes (Figure 4G) ----

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

## 2.1. Find all DE genes between the 2 conditions (part of Suppl. Table 7) ----
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

## 3.4. Vizualize GSEA results (Suppl. Figure 4G) ------------------------------
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
##### Here, we are not using the original LSC6 formula but Seurat's AddModuleScore, 
##### in order to account for the sparcity of the single cell data.
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

# 5. Cell type projections according to Zeng et al (2025) reference ============

##### Reference paper: https://doi.org/10.1158/2643-3230.BCD-24-0342
##### Reference Seurat object kindly provided by Andy Zeng

## 6.1. Load in the reference and look at it (plot not included in the paper) -----
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
               ~ FindTransferAnchors(reference = ref, query = .x, dims = 1:30, npcs = 30, reference.reduction = "pca"))

a1 <- FindTransferAnchors(reference = ref, query = seu_leukemic[[1]], dims = 1:30, npcs = 30, reference.reduction = "pca")
saveRDS(a1, paste0(wd, '377_anchrs_zeng_withRefRed_ptA.rds'))

a2 <- FindTransferAnchors(reference = ref, query = seu_leukemic[[2]], dims = 1:30, npcs = 30, reference.reduction = "pca")
saveRDS(a2, paste0(wd, '377_anchrs_zeng_withRefRed_ptB.rds'))

a3 <- FindTransferAnchors(reference = ref, query = seu_leukemic[[3]], dims = 1:30, npcs = 30, reference.reduction = "pca")
saveRDS(a3, paste0(wd, '377_anchrs_zeng_withRefRed_ptC.rds'))

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

##### also make Seurat objects in Andy's UMAP coordinates
seu_leukemic_umapAndy <- map2(anchors, seu_leukemic, ~ MapQuery(anchorset = .x, 
                         reference = ref, 
                         query = .y, 
                         refdata = list(Zeng_celltype = 'CellType_Broad'),
                         reference.reduction = "pca", 
                         reduction.model = "umap"))

pred_embeddings <- lapply(seu_leukemic_umapAndy, function(seuObj){
    pca <- Embeddings(seuObj[["ref.pca"]])
    umap <- Embeddings(seuObj[["ref.umap"]])
    out <- list('pca' = pca, 'umap' = umap)
    return(out)
  })

saveRDS(pred_embeddings, paste0(wd, '379_zengEmbeddings.rds'))

## 6.3. Visualize --------------------------------------------------------------

### 6.3.1. UMAP by predictions (Suppl. Figure 5B) 
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

### 6.3.2. UMAP by prediction confidence (Suppl. figure 6B)
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

### 6.3.3. Box of prediction confidence by condition (not included in the paper)
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

### 6.3.4. UMAP in Andy's coordinates, colored by cell density (not included in the paper)
ref_data <- Embeddings(subset(ref, downsample = 50000)[["umap"]]) %>%
  as.data.frame() 

query_data <- lapply(names(seu_leukemic_umapAndy), function(sample_name){
  seu_leukemic_umapAndy[[sample_name]]@meta.data$cell_label <- rownames(seu_leukemic_umapAndy[[sample_name]]@meta.data)
  df <- Embeddings(seu_leukemic_umapAndy[[sample_name]][['ref.umap']]) %>%
    as.data.frame() %>%
    mutate(patient = sample_name) %>%
    rename(UMAP_1 = refUMAP_1,
           UMAP_2 = refUMAP_2) %>%
    tibble::rownames_to_column(var = 'cell_label') %>%
    left_join(seu_leukemic_umapAndy[[sample_name]]@meta.data %>% 
                select(cell_label, predicted.Zeng_celltype.score),
              by = 'cell_label') %>%
    sample_n(3000)
  return(df)
})
names(query_data) <- names(seu_leukemic_umapAndy)

plotlist <- map2(query_data, names(query_data),
     ~ .x %>%
       ggplot(aes(x = UMAP_1, y = UMAP_2)) +
       geom_point(data = ref_data, color = '#E3E3E3', size = 0.05, alpha = 0.5) +
       ggpointdensity::geom_pointdensity(size = 0.2) +
       jcolors::scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1.75) +
       geom_density_2d(alpha = 0.4, color = 'black', h = 1.5, linewidth = 0.3) +
       theme_void() + 
       labs(subtitle = .y) +
       ggplot2::theme(strip.text.x = ggplot2::element_text(size = 18), 
                      legend.position = 'none'))

pdf(paste0(wd, '397_leukemicOnly_umapZeng_density.pdf'), height = 3.5, width = 10)
ggarrange(plotlist = plotlist, nrow = 1, ncol = 3)
dev.off()

### 6.3.5. UMAP in Andy's coordinates, colored by cell density, split by condition (Suppl. Figure 5C)
seu <- map(seu_leukemic_umapAndy, ~ SplitObject(.x, split.by = 'condition'))

query_data <- lapply(names(seu), function(MM_or_KD){
  sublist_new <- lapply(names(seu[[MM_or_KD]]), function(sample_name){
    seu[[MM_or_KD]][[sample_name]]@meta.data$cell_label <- rownames(seu[[MM_or_KD]][[sample_name]]@meta.data)
    df <- Embeddings(seu[[MM_or_KD]][[sample_name]][['ref.umap']]) %>%
      as.data.frame() %>%
      mutate(patient = sample_name) %>%
      rename(UMAP_1 = refUMAP_1,
             UMAP_2 = refUMAP_2) %>%
      tibble::rownames_to_column(var = 'cell_label') %>%
      left_join(seu[[MM_or_KD]][[sample_name]]@meta.data %>% 
                  select(cell_label, predicted.Zeng_celltype.score),
                by = 'cell_label') %>%
      sample_n(1500)
    return(df)
  })
  names(sublist_new) <- names(seu[[MM_or_KD]])
  return(sublist_new)
})
names(query_data) <- names(seu)

plotlist <- lapply(names(seu), function(patient_name){
  plts <- lapply(names(seu[[patient_name]]), function(MM_or_KD){
    plt <- query_data[[patient_name]][[MM_or_KD]] %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(data = ref_data, color = '#E3E3E3', size = 0.05, alpha = 0.5) +
      ggpointdensity::geom_pointdensity(size = 0.2) +
      jcolors::scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1.75) +
      geom_density_2d(alpha = 0.4, color = 'black', h = 1.5, linewidth = 0.3) +
      theme_void() + 
      labs(subtitle = MM_or_KD) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 18), 
                     legend.position = 'none')
    return(plt)
  })
  names(plts) <- names(seu[[patient_name]])
  return(plts)
})
names(plotlist) <- names(seu)

pdf(paste0(wd, '398_leukemicOnly_byCond_umapZeng_density.pdf'), height = 3.5, width = 6.5)

ggarrange(plotlist[['patientA']][['mismatch control']], 
          plotlist[['patientA']][['RUNX1::RUNX1T1 knockdown']], 
          nrow = 1, ncol = 2) +
  plot_annotation(title = 'Patient A')

ggarrange(plotlist[['patientB']][['mismatch control']], 
          plotlist[['patientB']][['RUNX1::RUNX1T1 knockdown']], 
          nrow = 1, ncol = 2) +
  plot_annotation(title = 'Patient B')

ggarrange(plotlist[['patientC']][['mismatch control']], 
          plotlist[['patientC']][['RUNX1::RUNX1T1 knockdown']], 
          nrow = 1, ncol = 2) +
  plot_annotation(title = 'Patient C (PDX)')

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
##### Available at https://doi.org/10.5281/zenodo.14578307, folder "scRNA"

# 98. Session info ==============================================================
sink(paste0(wd, '999_sessionInfo_3.txt'))
sessionInfo()
sink()