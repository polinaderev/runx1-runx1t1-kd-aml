library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(monocle3)

library(clusterProfiler)
library(msigdbr)

library(ggpubr)
library(patchwork)
library(viridis)
library(nVennR)

cbPalette <- c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "#999999")

cbPalette2 <- colorBlindness::paletteMartin
names(cbPalette2) <- NULL
cbPalette2 <- c(cbPalette2[3:9], cbPalette2[11:15], "#999999")

kdmm_palette <- c('RUNX1::RUNX1T1 knockdown' = '#f41626', 'mismatch control'= '#2538a5')

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
wd <- 'repos/runx1eto-kd-aml/scRNA/out/'

##### reproducibility
set.seed(42)

# 1. Load in data and prep it for integration ==================================

## 1.1. Load the Seurat objects with leukemic cells ----------------------------
seu_leukemic <- readRDS(paste0(wd, '300_seu_leukemic.rds'))

## 1.2. Replace orig.ident -----------------------------------------------------
origidents <- c('patientA', 'patientB', 'patientC')
names(origidents) <- names(seu_leukemic)

seu_leukemic <- map2(seu_leukemic, origidents, function(seuObj, identString){
  seuObj@meta.data$orig.ident <- as.character(seuObj@meta.data$orig.ident)
  seuObj@meta.data$orig.ident <- rep(identString, nrow(seuObj@meta.data))
  return(seuObj)
})

## 1.3. Downsample to the same number of cells per condition and per patient ----
targetCellNumber <- map(seu_leukemic, ~ table(.x@meta.data$condition)) %>%
  unlist() %>% 
  min()

seu_leukemic <- lapply(seu_leukemic, function(seuObj){
  Idents(seuObj) <- 'condition'
  return(seuObj)
})

seu_downsample <- map(seu_leukemic, ~ subset(.x, downsample = targetCellNumber))

# 2. Make one Seurat object out of the 3 patients' datasets ====================

## 2.1. Integrate --------------------------------------------------------------

features <- SelectIntegrationFeatures(object.list = seu_downsample, 
                                      nfeatures = 5000)
seu_downsample <- PrepSCTIntegration(object.list = seu_downsample, 
                                     anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seu_downsample, 
                                  normalization.method = "SCT", 
                                  anchor.features = features)

seu_integr <- IntegrateData(anchorset = anchors, 
                            normalization.method = "SCT")

## 2.2. Run dimensionality reduction on the integrated object ------------------

### 2.2.1 Run PCA
DefaultAssay(seu_integr) <- "integrated"
seu_integr <- RunPCA(seu_integr)

### 2.2.2 Vizualize PCA results

#### 2.2.2.1 Elbow plot (not included in the paper)
pdf(paste0(wd, "450_integr_elbow.pdf"))
ElbowPlot(seu_integr, ndims = 50) 
dev.off()

#### 2.2.2.2 PC heatmap (plot not included in the paper)
pdf(paste0(wd, "460_integr_dimheatmap.pdf"),
    width = 21,
    height = 21)
DimHeatmap(seu_integr, 
           dims = 1:50, 
           cells = 500, 
           balanced = TRUE)
dev.off()

#### 2.2.2.3. Scatter plot of the first 2 PCs, colored by condition (Supplementary Figure 4H)
pdf(paste0(wd, '470_integr_pca_byCond.pdf'), height = 3.5, width = 3)
DimPlot(seu_integr,
        reduction = 'pca',
        cols = kdmm_palette,
        group.by = "condition") +
  theme_void() +
  labs(title = '') +
  theme(legend.position = 'bottom')
dev.off()

##### Take 40 PCs for UMAP

### 2.2.3. Run UMAP

seu_integr <- RunUMAP(seu_integr, dims = 1:40) %>% 
  FindNeighbors(dims = 1:40)

saveRDS(seu_integr, paste0(wd, '475_leukemic_integrated.rds'))
##### seu_integr <- readRDS(paste0(wd, '475_leukemic_integrated.rds'))

### 2.2.4. Vizualize the UMAP

#### 2.2.4.1. UMAP by condition (Figure 4I)
pdf(paste0(wd, '480_integr_umap_byCond.pdf'), height = 3.5, width = 3)
DimPlot(seu_integr,
        reduction = 'umap',
        cols = kdmm_palette,
        group.by = "condition") +
  theme_void() +
  labs(title = '') +
  theme(legend.position = 'bottom')
dev.off()

#### 2.2.4.2. UMAP by patient (Figure 4J)
pdf(paste0(wd, '490_integr_umap_byPt.pdf'), height = 3.5, width = 3)
DimPlot(seu_integr,
        reduction = 'umap',
        cols = c(cbPalette[1], cbPalette[3], cbPalette[5]),
        group.by = "orig.ident",
        shuffle = TRUE) +
  theme_void() +
  labs(title = '') +
  theme(legend.position = 'bottom')
dev.off()

#### 2.2.4.3. UMAP by cell cycle stage (Supplementary Figure 6E)
pdf(paste0(wd, '500_integr_umap_byCellCycleStage.pdf'), height = 3.5, width = 3)
DimPlot(seu_integr,
        reduction = 'umap',
        cols = cbPalette,
        group.by = "Phase",
        shuffle = TRUE) +
  theme_void() +
  labs(title = '') +
  theme(legend.position = 'bottom')
dev.off()

# 3. Add the cell type predictions and embeddings (Zeng et al reference) and visualize ====

## 3.1. Add predictions and their scores ---------------------------------------

### 3.1.1. Load the cell type predictions
##### Produced in scRNA/03_leukemicCells_patientsSeparately.R
preds <- readRDS(paste0(wd, '378_zengPreds.rds'))

### 3.1.2. Fix the cell names
add <- c('1', '2', '3')

preds <- lapply(preds, function(sublist){
  sublist_new <- map2(sublist, add, function(vec, addition){
    names(vec) <- paste0(names(vec), '_', addition)
    return(vec)
  })
  vec_new <- c(sublist_new[[1]], sublist_new[[2]], sublist_new[[3]])
  return(vec_new)
})

### 3.1.3. Add the predictions and their scores to the Seurat object
seu_integr <- AddMetaData(seu_integr, metadata = preds[[1]], col.name = 'pred.Zeng.celltype') %>%
  AddMetaData(metadata = preds[[2]], col.name = 'pred.Zeng.score')

### 3.1.4. Visualize

#### 3.1.4.1. UMAP by cell type prediction (Figure 6A right)
pdf(paste0(wd, '506_integr_umap_byZengPred.pdf'), height = 3.5, width = 2.75)
DimPlot(seu_integr,
        reduction = 'umap',
        cols = zengPalette,
        group.by = "pred.Zeng.celltype",
        shuffle = TRUE) +
  theme_void() +
  labs(title = '') +
  theme(legend.position = 'bottom')
dev.off()

#### 3.1.4.2. UMAP by cell type prediction confidence (Figure 6B right)
pdf(paste0(wd, '514_integr_umap_byZengPred_score.pdf'), height = 3.5, width = 2.75)
FeaturePlot(seu_integr,
            features = 'pred.Zeng.score') +
  scale_color_viridis() +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom') +
  labs(title = NULL) 
dev.off()

#### 3.1.4.3. Box of cell type prediction confidence, by condition (Supplementary Figure 5A)
pdf(paste0(wd, '524_integr_box_byCond_ZengScore.pdf'), height = 2.5, width = 2)
ggplot(seu_integr@meta.data, 
       aes(x = condition, 
           y = pred.Zeng.score, 
           fill = condition, 
           colour = condition)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
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
                     label = "p.format")
dev.off()

## 3.2. Find markers for cell types and visualize ------------------------------

### 3.2.1. Find the markers
Idents(seu_integr) <- 'pred.Zeng.celltype'
DefaultAssay(seu_integr) <- 'RNA'

markers <- FindAllMarkers(seu_integr,
                          assay = 'RNA',
                          logfc.threshold = 0.1,
                          min.pct = 0.01,
                          only.pos = TRUE)
saveRDS(markers, paste0(wd, '530_zengCelltypes_markers.rds'))

### 3.2.2. Plot a heatmap with top markers -------------------------------------
markers <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.1 & p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1 | pct.2 > 0.1) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::slice_max(avg_log2FC, n = 10)

markers$cluster <- as.character(markers$cluster)
seu_integr@meta.data$pred.Zeng.celltype <- as.character(seu_integr@meta.data$pred.Zeng.celltype)

pdf(paste0(wd, '531_integr_heatmap_ZengCelltypes_markers.pdf'), height = 11, width = 15)
DoHeatmap(seu_integr %>% subset(subset = pred.Zeng.celltype %in% unique(markers$cluster)), 
          group.by = 'pred.Zeng.celltype',
          features = markers$gene,
          assay = "integrated",
          slot = "scale.data",
          group.bar = TRUE,
          group.colors = zengPalette) +
  scale_fill_viridis() +
  theme(legend.position = 'none')
dev.off()

## 3.3. Add the UMAP coordinates of my cells in Andy's UMAP coordinates and visualize ----

### 3.3.1. Load the embeddings (produced in scRNA/03_leukemicCells_patientsSeparately.R)
pred_embeddings <- readRDS(paste0(wd, '379_zengEmbeddings.rds'))

### 3.3.2. Fix the cell names
add <- c('1', '2', '3')

pred_embeddings <- map2(pred_embeddings, add, function(embedding_sublist, suffix){
  embedding_sublist_new <- lapply(embedding_sublist, function(embedding_df){
    rownames(embedding_df) <- paste0(rownames(embedding_df), '_', suffix)
    embedding_df <- embedding_df %>% as.data.frame() %>% rownames_to_column('cell_label')
    return(embedding_df)
  })
  return(embedding_sublist_new)
})

emb <- lapply(names(pred_embeddings[[1]]), function(pca_or_umap){
  df <- rbind(pred_embeddings[[1]][[pca_or_umap]],
              pred_embeddings[[2]][[pca_or_umap]],
              pred_embeddings[[3]][[pca_or_umap]])
  df <- df %>%
    filter(cell_label %in% colnames(seu_integr)) %>%
    column_to_rownames('cell_label')
  return(df)
})
names(emb) <- names(pred_embeddings[[1]])

seu_integr[['ref.pca']] <- CreateDimReducObject(as.matrix(emb[['pca']]))
seu_integr[['ref.umap']] <- CreateDimReducObject(as.matrix(emb[['umap']]))


# 4. Find the bone marrow cell type module scores ==============================

## 4.1 Find the necessary gene sets in MSigDB ----------------------------------
pathways <- msigdbr(species = 'Homo sapiens', category = 'C8')
pathways <- filter(pathways, grepl('^HAY_BONE_MARROW', gs_name))

mygenesets <- split(pathways, 
                    x = pathways$gene_symbol, 
                    f = pathways$gs_name)

## 4.2. Determine and plot the scores ------------------------------------------

### 4.2.1. Determine
DefaultAssay(seu_integr) <- 'RNA'
seu_integr <- AddModuleScore(seu_integr, features = mygenesets)
colnames(seu_integr@meta.data)[
  (length(colnames(seu_integr@meta.data))-length(mygenesets)+1):length(colnames(seu_integr@meta.data))
  ] <- 
  names(mygenesets)

### 4.2.2. Plot (Figure 7B)
pdf(paste0(wd, '542_integr_umap_byCond_celltypeScores_Hay_quantileColor.pdf'), height = 4, width = 6)
map(names(mygenesets),
    ~ FeaturePlot(seu_integr,
                  features = .x,
                  order = FALSE,
                  split.by = 'condition') &
      scale_colour_viridis(limits = c(quantile(seu_integr@meta.data[[.x]], probs = 0.05), 
                                      quantile(seu_integr@meta.data[[.x]], probs = 0.95)),
                           oob = scales::squish) &
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'bottom') &
      labs(title = NULL) &
      plot_annotation(subtitle = .x))
dev.off()

# 5. Plot some marker genes ====================================================
DefaultAssay(seu_integr) <- 'RNA'

## 5.1. UMAP not split (Figures 5C, S5E) ---------------------------------------
goi <- c('CD24', 'IL5RA', 'RETN', 'SERPINA1', 'RNASE2', 'CEBPE', 'GATA2')

pdf(paste0(wd, '545_integr_umap_genes_mono_granulo_eo_CD24_Cd125.pdf'), height = 3.5, width = length(goi)*2 + 2)
FeaturePlot(seu_integr,
            features = goi,
            order = TRUE,
            ncol = length(goi)) &
  scale_colour_viridis() &
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom')
dev.off()

## 5.2. UMAP split by condition (Figure 7A) ------------------------------------
goi <- 'CD34'
cd34_max <- max(
  as.data.frame(t(seu_integr$RNA@data))$CD34
)
cd34_min <- min(
  as.data.frame(t(seu_integr$RNA@data))$CD34
)

pdf(paste0(wd, '547_integr_umap_byCond_CD34.pdf'), height = 4, width = 6)
FeaturePlot(seu_integr,
                 features = goi,
                 order = TRUE,
                 split.by = 'condition',
                 combine = TRUE,
                 min.cutoff = cd34_min,
                 max.cutoff = cd34_max) &
  scale_colour_viridis(limits = c(cd34_min, cd34_max)) &
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom')
dev.off()

# 6. Do the pseudotime analysis ================================================

## 6.1. Cluster and make the graph ---------------------------------------------
cds <- as.cell_data_set(seu_integr) 

cds <- cluster_cells(cds, k = 250, cluster_method = 'louvain') %>%
  learn_graph()

## 6.2. Root the graph ---------------------------------------------------------

### 6.2.1. Which cells are the 1% most immature cells?
quant99 <- quantile(seu_integr@meta.data$HAY_BONE_MARROW_CD34_POS_HSC, probs = 0.99)
immature_cells <- rownames(seu_integr@meta.data %>% 
                             dplyr::filter(HAY_BONE_MARROW_CD34_POS_HSC > quant99))

### 6.2.2. Pick the graph node that is the most heavily occupied by immature cells
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])

root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[immature_cells,]))))]

### 6.2.3. Root the graph
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

## 6.3. Plot the pseudotime (Figure 4K) ----------------------------------------
pdf(paste0(wd, '550_integr_umap_trajectory_byMonoclePseudotime.pdf'), 
    width = 3, height = 3.5)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_cell_groups = FALSE,
           show_trajectory_graph = TRUE) +
  scale_colour_viridis() +
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom')
dev.off()

# 7. Low resolution clustering analysis (superclusters) ========================

## 7.1. Cluster and split by condition and by supercluster ---------------------
DefaultAssay(seu_integr) <- 'integrated'

seu_integr <- FindClusters(seu_integr, resolution = 0.025)

seu_cond <- SplitObject(seu_integr, split.by = 'condition')
seu_superclust <- SplitObject(seu_integr, split.by = 'integrated_snn_res.0.025')

## 7.2. Plot the superclusters (Figure 6C) -------------------------------------
pdf(paste0(wd, '560_integr_umap_2clust.pdf'), width = 3, height = 3.5)
map2(seu_cond, names(seu_cond),
     ~ DimPlot(.x,
               reduction = 'umap',
               group.by = 'integrated_snn_res.0.025',
               cols = c(cbPalette[6], cbPalette[7]),
               label = TRUE) +
       theme_void() +
       labs(title = NULL, subtitle = .y) +
       theme(legend.position = 'none'))
dev.off()

## 7.3. Boxplot by cell type prediction confidence by supercluster (Figure 6E) -----
pdf(paste0(wd, '570_integr_box_by2Clust_ZengScore.pdf'), height = 2.5, width = 2.5)
map2(seu_superclust, names(seu_superclust),
     ~ ggplot(.x@meta.data, 
              aes(x = condition, 
                  y = pred.Zeng.score, 
                  fill = condition, 
                  colour = condition)) +
       geom_violin() +
       geom_boxplot(width = 0.1) +
       theme_minimal() +
       scale_fill_manual(values = kdmm_palette) +
       scale_color_manual(values = c('black', 'black')) +
       labs(y = 'Zeng score',
            x = 'cluster',
            subtitle = .y) +
       theme(legend.position = 'none',
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "black"),
             axis.ticks = element_line(colour = "black"),
             axis.text.x = element_text(colour = "black"),
             axis.text.y = element_text(colour = "black")) +
       stat_compare_means(comparisons = list(c(1, 2)), 
                          method = 'wilcox.test',
                          label = "p.format"))
dev.off()

## 7.4. Proportions of phenotypes in supercluster I (Figure 6D left) ----------------
seu_siRE_superclust <- SplitObject(seu_cond[['RUNX1::RUNX1T1 knockdown']], split.by = 'integrated_snn_res.0.025')

celltype_counts <- dplyr::count(seu_siRE_superclust[['0']]@meta.data, 
                                pred.Zeng.celltype) %>%
  mutate(condition = 'knockdown')

celltype_counts$pred.Zeng.celltype <- factor(celltype_counts$pred.Zeng.celltype,
                                                 levels = unique(celltype_counts$pred.Zeng.celltype))

pdf(paste0(wd, '580_integr_kd_area_2clust_clust0_ZengPreds.pdf'), width = 3, height = 2)
ggplot(celltype_counts,
       aes(fill = pred.Zeng.celltype,
           y = n,
           x = condition)) +
  geom_col(position = "fill") +
  labs(x = 'condition',
       y = 'proportion of phenotype') +
  scale_fill_manual(values = zengPalette) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
dev.off()

## 7.5. Proportions of phenotypes in supercluster II (Figure 6D right) ---------
celltype_counts <- dplyr::count(seu_siRE_superclust[['1']]@meta.data, 
                                pred.Zeng.celltype) %>%
  mutate(condition = 'knockdown')

celltype_counts$pred.Zeng.celltype <- factor(celltype_counts$pred.Zeng.celltype,
                                             levels = unique(celltype_counts$pred.Zeng.celltype))

pdf(paste0(wd, '588_integr_kd_area_2clust_clust1_ZengPreds.pdf'), width = 3, height = 2)
ggplot(celltype_counts,
       aes(fill = pred.Zeng.celltype,
           y = n,
           x = condition)) +
  geom_col(position = "fill") +
  labs(x = 'condition',
       y = 'proportion of phenotype') +
  scale_fill_manual(values = zengPalette) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
dev.off()

## 7.7. Cell type module scores (calculated earlier, now plotted for KD only; Supplementary Figure 5B) -----
pdf(paste0(wd, '600_integr_kd_umap_HayModuleScores_quantileColor.pdf'), height = 4, width = 3)
map(names(mygenesets),
    ~ FeaturePlot(seu_cond[['RUNX1::RUNX1T1 knockdown']],
                  features = .x,
                  order = FALSE,
                  max.cutoff = 'q95',
                  min.cutoff = 'q5') +
      scale_colour_viridis() +
      theme(axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = 'bottom') +
      labs(title = NULL) +
      plot_annotation(subtitle = .x))
dev.off()

## 7.8. Find markers of the clusters and compare them with markers of conditions ----

### 7.8.1. Markers of the clusters (Supplementary Table 8)
Idents(seu_integr) <- 'integrated_snn_res.0.025'

markers_clust <- FindMarkers(seu_integr,
                       ident.1 = '0',
                       assay = 'RNA',
                       logfc.threshold = 0,
                       min.pct = 0)

saveRDS(markers_clust, paste0(wd, '601_integr_markers_clust0_vs_clust1.rds'))
##### markers_clust <- readRDS(paste0(wd, '601_integr_markers_clust0_vs_clust1.rds'))
write.csv(markers_clust, 
          paste0(wd, '601_integr_markers_clust0_vs_clust1.csv'),
          row.names = TRUE)

### 7.8.2. Markers of the conditions (Supplementary Table 7)
Idents(seu_integr) <- 'condition'

markers_cond <- FindMarkers(seu_integr,
                       ident.1 = 'RUNX1::RUNX1T1 knockdown',
                       assay = 'RNA',
                       logfc.threshold = 0,
                       min.pct = 0)

saveRDS(markers_cond, paste0(wd, '602_integr_markers_REKD_vs_MM.rds'))
##### markers_cond <- readRDS(paste0(wd, '602_integr_markers_REKD_vs_MM.rds'))
write.csv(markers_cond, 
          paste0(wd, '602_integr_markers_REKD_vs_MM.csv'),
          row.names = TRUE)

### 7.8.3. Make a list with the 2 tables of markers
markers <- list(
  'clust' = markers_clust,
  'cond' = markers_cond
)

### 7.8.4. nVenn diagram of markers of the clusters and conditions (Figure 6F)
de_pos <- lapply(markers, function(tib){
  tib$gene <- rownames(tib)
  tib_pos <- dplyr::filter(tib, avg_log2FC > 0 & p_val_adj < 0.001)
  genes_pos <- tib_pos$gene
  return(genes_pos)
})

de_neg <- lapply(markers, function(tib){
  tib$gene <- rownames(tib)
  tib_neg <- dplyr::filter(tib, avg_log2FC < 0 & p_val_adj < 0.001)
  genes_neg <- tib_neg$gene
  return(genes_neg)
})

markers_venn <- list(
  'Up in cluster 0' = de_pos[['clust']],
  'Up in cluster 1' = de_neg[['clust']],
  'Up in siRE' = de_pos[['cond']],
  'Up in siMM' = de_neg[['cond']]
)

plotVenn(markers_venn, 
         showPlot = TRUE, 
         systemShow = FALSE, 
         outFile = paste0(wd, '603_integr_nvenn_clust_cond.svg'),
         setColors = c(cbPalette[6], cbPalette[7], kdmm_palette[1], kdmm_palette[2])
         )

### 7.8.5. Lists of genes in each area of the Nvenn diagram (not included in the paper but you can re-make it by overlapping Suppl. Tables 7 and 8)
result <- c()

groups <- names(upregulated_genes)
for (i in seq_along(groups)) {
  for (j in (i + 1):length(groups)) {
    group1 <- groups[i]
    group2 <- groups[j]
    overlap <- intersect(upregulated_genes[[group1]], upregulated_genes[[group2]])
    if (length(overlap) > 0) {
      overlap_genes <- paste(overlap, collapse = ", ")
      result <- c(result, paste0(group1, " & ", group2, ": ", overlap_genes))
    }
  }
}

exclusive_conditions <- list(
  "Up in cluster 0 & NOT Up in siMM & NOT Up in siRE" = setdiff(
    upregulated_genes[["Up in cluster 0"]],
    c(upregulated_genes[["Up in siMM"]], upregulated_genes[["Up in siRE"]])
  ),
  "Up in cluster 1 & NOT Up in siMM & NOT Up in siRE" = setdiff(
    upregulated_genes[["Up in cluster 1"]],
    c(upregulated_genes[["Up in siMM"]], upregulated_genes[["Up in siRE"]])
  ),
  "Up in siMM & NOT Up in cluster 0 & NOT Up in cluster 1" = setdiff(
    upregulated_genes[["Up in siMM"]],
    c(upregulated_genes[["Up in cluster 0"]], upregulated_genes[["Up in cluster 1"]])
  ),
  "Up in siRE & NOT Up in cluster 0 & NOT Up in cluster 1" = setdiff(
    upregulated_genes[["Up in siRE"]],
    c(upregulated_genes[["Up in cluster 0"]], upregulated_genes[["Up in cluster 1"]])
  )
)

for (condition in names(exclusive_conditions)) {
  genes <- exclusive_conditions[[condition]]
  if (length(genes) > 0) {
    result <- c(result, paste0(condition, ": ", paste(genes, collapse = ", ")))
  }
}

writeLines(result, paste0(wd, "604_geneListsFor_nVenn.txt"))

## 7.9. UMAP split by condition, colored by patient and with nr of cells per patient in each supercluster (Suppl. Fig. 5) ------------
seu_cond <- SplitObject(seu_integr, split.by = 'condition')

pdf(paste0(wd, '604_integr_umap_2clust_byCond_byPt.pdf'), width = 3, height = 3.5)
map2(seu_cond, names(seu_cond),
     ~ DimPlot(.x,
               reduction = 'umap',
               group.by = 'orig.ident',
               cols = c(cbPalette[1], cbPalette[3], cbPalette[5]),
               shuffle = TRUE) +
       theme_void() +
       labs(title = NULL, subtitle = .y) +
       theme(legend.position = 'bottom'))
dev.off()

map(seu_cond, ~ table(.x@meta.data$integrated_snn_res.0.025, .x@meta.data$orig.ident))

# $`RUNX1::RUNX1T1 knockdown`
# 
# patientA patientB patientC
# 0      997      745     1341
# 1      552      804      208
# 
# $`mismatch control`
# 
# patientA patientB patientC
# 0      839      941     1362
# 1      710      608      187


# 8. High resolution clustering analysis =======================================

## 8.1. Cluster ----------------------------------------------------------------
DefaultAssay(seu_integr) <- 'integrated'

seu_integr <- FindClusters(seu_integr, resolution = 0.57)

seu <- SplitObject(seu_integr, split.by = 'condition')

## 8.2. Plot the clusters (Figure 5A) ------------------------------------------
pdf(paste0(wd, '605_integr_umap_hiResClust.pdf'), width = 3, height = 3.5)
map2(seu, names(seu),
     ~ DimPlot(.x,
               reduction = 'umap',
               group.by = 'integrated_snn_res.0.57',
               cols = cbPalette2,
               label = TRUE) +
       theme_void() +
       labs(title = NULL, subtitle = .y) +
       theme(legend.position = 'none'))
dev.off()

## 8.3. Find top markers for each cluster (not included in the paper) ---------
Idents(seu_integr) <- 'integrated_snn_res.0.57'
markers <- FindAllMarkers(seu_integr, 
                          assay = "RNA",
                          slot = "data",
                          only.pos = FALSE, 
                          min.pct = 0.1, 
                          logfc.threshold = 0.1) %>%
  dplyr::group_by(cluster) 

saveRDS(markers, paste0(wd, '610_integr_highRes_markers.rds'))
write_csv(markers, paste0(wd, '610_integr_highRes_markers.csv'))

## 8.4. Plot a heatmap with top markers (Figure 5B) ----------------------------
markers <- dplyr::filter(markers, avg_log2FC > 0.1 & p_val_adj < 0.05) %>%
  dplyr::filter(pct.1 > 0.1 | pct.2 > 0.1) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::slice_max(avg_log2FC, n = 10)

pdf(paste0(wd, '621_integr_heatmap_hiResClust_markers.pdf'), height = 8, width = 5)
DoHeatmap(seu_integr, 
          group.by = 'integrated_snn_res.0.57',
          features = markers$gene,
          assay = "integrated",
          slot = "scale.data",
          group.bar = TRUE,
          group.colors = cbPalette2) +
  scale_fill_viridis()
dev.off()

# 10. GSEA =====================================================================

## 10.1. Get the gene sets ------------------------------------------------------
pathways <- msigdbr(species = 'Homo sapiens', category = 'C8')
goi <- filter(pathways, 
              gs_name == 'ZHENG_CORD_BLOOD_C6_HSC_MULTIPOTENT_PROGENITOR'|gs_name == 'HAY_BONE_MARROW_CD34_POS_HSC') %>%
  select(gs_name, gene_symbol)

## 10.2. Make ranks ------------------------------------------------------------
markers_cond_filt <- markers_cond %>%
  filter(pct.1 + pct.2 > 0.01) ##### filter out genes whose expression is low
ranks <- markers_cond_filt$avg_log2FC
names(ranks) <- rownames(markers_cond_filt)
ranks <- ranks[!is.na(ranks)] %>% 
  sort(decreasing = TRUE) ##### sorting is required for ClusterProfiler

## 10.3. Run GSEA --------------------------------------------------------------
gsea_res <- GSEA(
  geneList = ranks,
  maxGSSize = 2000,
  pvalueCutoff = 1,
  eps = 0,
  pAdjustMethod = 'BH',
  seed = TRUE,
  TERM2GENE = goi
)

## 10.4. Vizualize GSEA results (Figure 7C, Supplementary Figure 6A) -----------
pdf(paste0(wd, '640_integr_gseaPlots_OlafsRanks_Zheng_Hay_HSC.pdf'), width = 10, height = 10)
plotlist <- lapply(paste0(unique(goi$gs_name)), function(geneset_name){
  p <- enrichplot::gseaplot(
    gsea_res,
    geneSetID = geneset_name,
    title = paste0('Integrated ptABC', ', ', geneset_name),
    color.line = '#21908c') %>%
    as.patchwork()
  p <- p + plot_annotation(subtitle = paste0(
    'NES = ', as.character(signif(gsea_res@result$NES[gsea_res@result$ID == geneset_name], digits = 3)),
    ', padj = ', as.character(signif(gsea_res@result$p.adjust[gsea_res@result$ID == geneset_name], digits = 3))
  ))
  return(p)
})
print(plotlist[[1]])
print(plotlist[[2]])
dev.off()

# 11. LSC 6 score (Figure 7E) ==============================================================
LSC_signature <- c("DNMT3B", "CD34", "ADGRG1", "SOCS2", "SPINK2", "FAM30A")

seu_integr <- AddModuleScore(seu_integr, 
                             features = list(LSC_signature),
                             assay = 'RNA',
                             name = 'LSC6')

seu_integr@meta.data <- rename(seu_integr@meta.data, LSC6 = LSC61)

pdf(paste0(wd, '651_integr_umap_byCond_pediatricLSC6.pdf'), height = 3.5, width = 6)
FeaturePlot(seu_integr,
            features = "LSC6",
            order = TRUE,
            split.by = 'condition') &
  scale_colour_viridis(limits = c(min(seu_integr@meta.data$LSC6), 
                                  max(seu_integr@meta.data$LSC6))) &
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'bottom') &
  labs(title = NULL)

dev.off()

# 13. Make supplementary tables 7 and 8 ========================================

## 13.1. Supplementary table 7: markers RE vs MM -------------------------------
markers <- readRDS(paste0(wd, '328_leukemic_allMarkers_KDvsMM.rds'))
int <- readRDS(paste0(wd, '602_integr_markers_REKD_vs_MM.rds'))

markers[[4]] <- int
names(markers) <- c('ptA', 'ptB', 'ptC', 'ptABCintegrated')
markers <- lapply(names(markers), function(name){
  df <- markers[[name]]
  colnames(df)[3] <- 'siRE_proportionPositiveCells'
  colnames(df)[4] <- 'siMM_proportionPositiveCells'
  colnames(df) <- paste0(name, '_', colnames(df))
  df$gene <- rownames(df)
  return(df)
})

markers_df <- markers[[1]] %>%
  left_join(markers[[2]], by = 'gene') %>%
  left_join(markers[[3]], by = 'gene') %>%
  left_join(markers[[4]], by = 'gene')

write_csv(markers_df, paste0(wd, '680_allMarkers_allPatients_RE_vs_MM.csv'))

## 13.2. Supplementary table 8: markers supercluster I vs II -------------------
full <- readRDS(paste0(wd, '601_integr_markers_clust0_vs_clust1.rds'))
re <- readRDS(paste0(wd, '670_integr_REonly_markers_clust0_vs_clust1.rds'))

markers <- list('ptABC_siRE_and_siMM' = full, 'ptABC_siREonly' = re)
markers <- lapply(names(markers), function(name){
  df <- markers[[name]]
  colnames(df)[3] <- 'superclustI_proportionPositiveCells'
  colnames(df)[4] <- 'superclustII_proportionPositiveCells'
  colnames(df) <- paste0(name, '_', colnames(df))
  df$gene <- rownames(df)
  return(df)
})

markers_df <- markers[[1]] %>%
  left_join(markers[[2]], by = 'gene')

write_csv(markers_df, paste0(wd, '690_allMarkers_superclustI_vs_superclustII.csv'))


# 98. Make a Seurat object to upload to Zenodo =================================
seu_save <- seu_integr
seu_save[['HTO']] <- NULL
seu_save[['SCT']] <- NULL
seu_save@meta.data <- seu_save@meta.data %>%
  select(orig.ident, Phase, condition, pred.Zeng.celltype, pred.Zeng.score,
         pred.vanGalen.celltype, pred.vanGalen.score, matches('^HAY_BONE_MARROW'),
         matches('^integrated_snn_res'))
saveRDS(seu_save, paste0(wd, 'seu_ptABCintegrated_leukemicCells.rds'))

# 99. Session info =============================================================
sink(paste0(wd, '999_sessionInfo_4.txt'))
sessionInfo()
sink()