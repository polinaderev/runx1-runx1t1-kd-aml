library(Seurat)
library(tidyverse)

library(patchwork)
library(viridis)
library(ggpubr)
library(Nebulosa)

##### The color-blind-friendly palette:
cbPalette <- c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "#999999")

##### Another color-blind-friendly palette:
cbPalette2 <- colorBlindness::paletteMartin
names(cbPalette2) <- NULL
cbPalette2 <- c(cbPalette2[3:9], cbPalette2[11:15], "#999999")

##### Directory for outputs:
wd <- 'repos/runx1eto-kd-aml/scRNA/out/'

##### reproducibility
set.seed(42)

# 1. Load in data ==============================================================

## 1.1. UMI and HTO count matrix (from cellranger count output) ----------------
in_dir <- 'EX113_LS/deep/01_cellranger/'
patient_name <- c('A', 'B', 'C')

seu <- lapply(patient_name, function(name){
  dat <- Read10X(data.dir = paste0(in_dir, 'pt', name, '_deep/outs/filtered_feature_bc_matrix'))
  umis <- dat$`Gene Expression`
  htos <- dat$`Antibody Capture`
  joint.bcs <- intersect(colnames(umis), colnames(htos))
  umis <- umis[, joint.bcs]
  htos <- as.matrix(htos[, joint.bcs])
  seuObj <- CreateSeuratObject(counts = umis)
  lst <- list('seuObj' = seuObj, 'htos' = htos)
  return(lst)
})
names(seu) <- paste0('patient', patient_name)

### 1.1.1. Fix HTO row names for patient C
rownames(seu[['patientC']][['htos']])[rownames(seu[['patientC']][['htos']]) == 'hto10-mm'] <- 'hto-mm' ##### mm = mismatch control
rownames(seu[['patientC']][['htos']])[rownames(seu[['patientC']][['htos']]) == 'hto9-re'] <- 'hto-kd' ##### kd = RUNX1/ETO knockdown

## 1.2. SoupOrCell donor info (from souporcell output) -------------------------
in_dir <- 'EX113_LS/deep/02_souporcell/'
patient_name <- c('A', 'B', 'C')

donors <- lapply(patient_name, function(name){
  tib <- read_tsv(paste0(in_dir, 'pt', name, '/clusters.tsv'))
  return(tib)
})
names(donors) <- names(seu)

# 2. Demultiplex ===============================================================

## 2.1. Normalize --------------------------------------------------------------
seu <- lapply(seu, function(sublist){
  new_seuObj <- NormalizeData(sublist$seuObj) %>%
    FindVariableFeatures(selection.method = "mean.var.plot") %>%
    ScaleData(features = VariableFeatures(sublist$seuObj))
  new_list <- list('seuObj' = new_seuObj, 'htos' = sublist$htos)
  return(new_list)
})

## 2.2 Add HTO data as a new assay independent from RNA ------------------------
##### Patients A and B have 2 hashtags per condition, patient C only has one.
seu_ptC <- seu$patientC
seu_ptAB <- seu[c(1,2)]

seu_ptAB <- lapply(seu_ptAB, function(sublist){
  hto_mm <- sublist$htos[1,] + sublist$htos[2,] ##### mm = mismatch control
  hto_kd <- sublist$htos[3,] + sublist$htos[4,] ##### kd = RUNX1::RUNX1T1 knockdown
  htos_merged <- rbind(hto_mm, hto_kd)
  sublist$seuObj[['HTO']] <- CreateAssayObject(counts = htos_merged)
  return(sublist$seuObj)
})

seu_ptC$seuObj[['HTO']] <- CreateAssayObject(counts = seu_ptC$htos)
seu_ptC <- seu_ptC$seuObj

seu <- append(seu_ptAB, seu_ptC)
names(seu)[3] <- 'patientC'

## 2.3. Normalize HTO data, here centered log-ratio transformation is used -----
seu <- map(seu, 
           ~ NormalizeData(.x, 
                           assay = "HTO", 
                           normalization.method = "CLR", 
                           margin = 1))

## 2.4. Demultiplex cells based on HTO enrichment ------------------------------
seu <- map(seu,
           ~ HTODemux(.x, 
                      assay = "HTO", 
                      positive.quantile = 0.95))

##### look at the results
lapply(seu, function(seuObj){
  table(seuObj$HTO_classification.global)
})

##### try another positive quantile
seu <- map(seu,
           ~ HTODemux(.x, 
                      assay = "HTO", 
                      positive.quantile = 0.99))

##### look at the results
lapply(seu, function(seuObj){
  table(seuObj$HTO_classification.global)
})

##### I prefer the 0.95 positive quantile.
seu <- map(seu,
           ~ HTODemux(.x, 
                      assay = "HTO", 
                      positive.quantile = 0.95))

# 2.5. Demultiplex based on SoupOrCell assignment ------------------------------
donors <- lapply(donors, function(df){
  newnames <- paste0('soup_', colnames(df))
  df <- rename_at(df, vars(colnames(df)), ~ newnames)
  rownames(df) <- df$soup_barcode
  return(df)
})

seu <- map2(seu, donors,
            function(obj, new_metadata){
              obj@meta.data$soup_barcode <- rownames(obj@meta.data)
              obj@meta.data <- left_join(obj@meta.data, new_metadata, by = 'soup_barcode')
              rownames(obj@meta.data) <- obj@meta.data$soup_barcode 
              obj@meta.data <- dplyr::select(obj@meta.data, -soup_barcode)
              return(obj)
            })

## 2.6. Vizualize --------------------------------------------------------------
seu <- lapply(seu, function(seuObj){
  Idents(seuObj) <- 'HTO_maxID'
  return(seuObj)
})

### 2.6.1. 010 HTO ridgeplot (plot not included in the paper)
pdf(paste0(wd, '010_HTO_ridgeplot.pdf'), height = 3)
map2(seu, names(seu),
     ~ RidgePlot(.x, 
          assay = 'HTO', 
          features = rownames(.x[['HTO']]),
          cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

### 2.6.2. 020 Pairs of HTO signals to confirm mutual exclusivity in singlets (plot not included in the paper)
pdf(paste0(wd, '020_HTO_scatter.pdf'))
map2(seu, names(seu),
     ~ FeatureScatter(.x, 
                      feature1 = "hto-mm", 
                      feature2 = "hto-kd",
                      cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

### 2.6.3. 030 Compare number of UMIs for singlets, doublets and negative cells (plot not included in the paper) 
seu <- lapply(seu, function(seuObj){
  Idents(seuObj) <- 'HTO_classification.global'
  return(seuObj)
})

pdf(paste0(wd, "030_HTO_violin.pdf"))
map2(seu, names(seu),
     ~ VlnPlot(.x, 
        features = "nCount_RNA", 
        pt.size = 0.1, 
        log = TRUE,
        cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

### 2.6.4. 040 HTO heatmap (plot not included in the paper)
pdf(paste0(wd, "040_HTO_heatmap.pdf"))
map2(seu, names(seu),
     ~ HTOHeatmap(.x, 
           assay = "HTO") +
       scale_fill_viridis() + 
    plot_annotation(title = .y))
dev.off()

### 2.6.5. 050 Preliminary RNA UMAP colored by HTO assignment (plot not included in the paper)
seu <- lapply(seu, function(seuObj){
  DefaultAssay(seuObj) <- 'RNA'
  return(seuObj)
})

seu <- map(seu, 
           ~ SCTransform(.x) %>%
             RunPCA() %>%
             RunUMAP(dims = 1:40) %>%
             FindNeighbors(dims = 1:40))

pdf(paste0(wd, "050_rna_umap_byHTO.pdf"))
map2(seu, names(seu),
     ~ DimPlot(.x,
               group.by = "hash.ID",
               cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

### 2.6.6. 060 RNA UMAP colored by souporcell assignment (doublet, singlet or negative) (plot not included in the paper)
pdf(paste0(wd, "060_rna_umap_bySoup.pdf"))
map2(seu, names(seu),
     ~ DimPlot(.x,
        group.by = "soup_status",
        cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

### 2.6.7. 070 RNA UMAP colored by souporcell assignment (which donor) (plot not included in the paper)
pdf(paste0(wd, "070_rna_umap_bySoupAssignment.pdf"))
map2(seu, names(seu),
     ~ DimPlot(.x,
               group.by = "soup_assignment",
               cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

## 2.7. Get rid of doublets and unassigned droplets ----------------------------
seu <- lapply(seu, function(obj){
  obj_filt <- subset(obj, soup_assignment == '0' | soup_assignment == '1') %>%
    subset(soup_status == 'singlet' & HTO_classification.global == 'Singlet') 
  return(obj_filt)
})

## 2.8. Vizualize the filtered object ------------------------------------------

### 2.8.1. 080 RNA UMAP colored by HTO assignment (plot not included in the paper)
pdf(paste0(wd, "080_filt_rna_umap_byHTO.pdf"))
map2(seu, names(seu),
     ~ DimPlot(.x,
               group.by = "hash.ID",
               cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

### 2.8.2. 090 RNA UMAP colored by souporcell assignment (doublet, singlet or negative) (plot not included in the paper)
pdf(paste0(wd, "090_filt_rna_umap_bySoup.pdf"))
map2(seu, names(seu),
     ~ DimPlot(.x,
               group.by = "soup_status",
               cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

### 2.8.3. 100 RNA UMAP colored by souporcell assignment (which donor) (plot not included in the paper)
pdf(paste0(wd, "100_filt_rna_umap_bySoupAssignment.pdf"))
map2(seu, names(seu),
     ~ DimPlot(.x,
               group.by = "soup_assignment",
               cols = cbPalette) + 
       plot_annotation(title = .y))
dev.off()

# 3. Filter for cells/droplets of low quality ==================================

## 3.1 Add percentage of mitochondrial transcripts to the metadata -------------

seu <- lapply(seu, function(seuObj){
  seuObj[['percent.mt']] <- PercentageFeatureSet(seuObj, pattern = '^MT-')
  return(seuObj)
})

## 3.2 Vizualize the QC metrics ------------------------------------------------

### 3.2.1. Violin plot: features, counts, mitochondrial transcript percentage (plot not included in the paper)
seu <- lapply(seu, function(seuObj){
  Idents(seuObj) <- 'orig.ident'
  return(seuObj)
})

pdf(paste0(wd, '110_vln_qcmetr.pdf'))
map2(seu, names(seu),
     ~ VlnPlot(.x, 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 3,
               cols = cbPalette) +
       plot_annotation(title = .y)
)
dev.off()

### 3.2.2. Violin: mitochondrial percentage on the range from 0 to 20% (plot not included in the paper)
pdf(paste0(wd, '120_vln_percentmt0-20.pdf'), width = 2)
map2(seu, names(seu),
     ~ VlnPlot(.x, features = "percent.mt") +
       ylim(0, 20) +
       plot_annotation(title = .y) +
       theme(legend.position='none'))
dev.off()

### 3.2.3. Scatter: "percent.mt" and "nFeature_RNA" (plot not included in the paper)
pdf(paste0(wd, '130_scatter_nfeature_percentmt.pdf'))
map2(seu, names(seu),
     ~FeatureScatter(.x, 
                     feature1 = "nFeature_RNA", 
                     feature2 = "percent.mt") +
       plot_annotation(title = .y)
)
dev.off()

### 3.2.4. UMAP by mitochondrial percentage (plot not included in the paper)
pdf(paste0(wd, '140_umap_percentmt.pdf'))
map2(seu, names(seu),
     ~FeaturePlot(.x, features = "percent.mt") + 
       scale_color_viridis() +
       plot_annotation(title = .y)
)
dev.off()

### 3.2.5. UMAP by nfeature (plot not included in the paper)
pdf(paste0(wd, '150_umap_nfeature.pdf'))
map2(seu, names(seu),
     ~FeaturePlot(.x, features = "nFeature_RNA") + 
       scale_color_viridis() +
       plot_annotation(title = .y)
)
dev.off()

### 3.2.6. UMAP by ncount (plot not included in the paper)
pdf(paste0(wd, '160_umap_ncount.pdf'))
map2(seu, names(seu),
     ~FeaturePlot(.x, features = "nCount_RNA") + 
       scale_color_viridis() +
       plot_annotation(title = .y)
)
dev.off()

### 3.2.7. Find clusters to determine if cells with low mitochondrial percentage belong to the same cluster
seu <- lapply(seu, function(seuObj){
  DefaultAssay(seuObj) <- 'SCT'
  return(seuObj)
})

seu <- map(seu,
           ~ FindNeighbors(.x) %>% FindClusters())

### 3.2.8. UMAP by cluster (plot not included in the paper)
pdf(paste0(wd, '170_umap_cluster.pdf'))
map2(seu, names(seu),
     ~ DimPlot(.x) +
       plot_annotation(title = .y))
dev.off()

### 3.2.9 Violin of nfeatures on the range from 0 to 2000 (plot not included in the paper)
seu <- lapply(seu, function(seuObj){
  Idents(seuObj) <- 'HTO_classification.global'
  return(seuObj)
})

pdf(paste0(wd, '180_vln_nfeature_0-2000.pdf'))
map2(seu, names(seu),
     ~ VlnPlot(.x, features = "nFeature_RNA") +
       ylim(0, 2000) +
       plot_annotation(title = .y) +
       theme(legend.position = 'none'))
dev.off()

##### Conclusion from the QC plots: 7% mitochondrial transcripts and 750 features should be the cutoff.

## 3.3. Filter -----------------------------------------------------------------

seu <- map(seu, ~ subset(.x, 
                         percent.mt < 7 & nFeature_RNA > 750))

# 4. Dimensionality reduction ==================================================

## 4.1. Run SCTransform, regress for difference between G2M and S. Run PCA -----
seu <- lapply(seu, function(seuObj){
  DefaultAssay(seuObj) <- 'RNA'
  seuObj <-  CellCycleScoring(seuObj, 
                              s.features = cc.genes$s.genes, 
                              g2m.features = cc.genes$g2m.genes)
  seuObj$CC.Difference <- seuObj$S.Score - seuObj$G2M.Score
  seuObj <- SCTransform(seuObj, vars.to.regress = 'CC.Difference') %>% RunPCA()
  return(seuObj)
})

## 4.2. Determine how many PCs to take ------------------------------------------

### 4.2.1. Elbow plot (plot not included in the paper)
pdf(paste0(wd, '190_elbow.pdf'))
map2(seu, names(seu),
     ~ ElbowPlot(.x,
                ndims = 50) +
       plot_annotation(title = .y))
dev.off()

### 4.2.2. PC heatmap (plot not included in the paper)
pdf(paste0(wd, '200_dimheatmap.pdf'),
    width = 32,
    height = 21)
map2(seu, names(seu),
    ~ DimHeatmap(.x, 
                 dims = 1:50, 
                 cells = 500, 
                 balanced = TRUE,
                 fast = FALSE,
                 ncol = 5) +
      plot_annotation(title = .y) &
      scale_fill_viridis() 
)
dev.off()

###### Patient A: take 35 PCs
###### Patient B: take 29 PCs
###### Patient C: take 19 PCs

## 4.3. Find neighbors and run UMAP with the abovementioned amounts of PCs -----

seu[['patientA']] <- FindNeighbors(seu[['patientA']], dims = 1:35) %>% RunUMAP(dims = 1:35)
seu[['patientB']] <- FindNeighbors(seu[['patientB']], dims = 1:29) %>% RunUMAP(dims = 1:29)
seu[['patientC']] <- FindNeighbors(seu[['patientC']], dims = 1:19) %>% RunUMAP(dims = 1:19)

## 4.4. Vizualize results of dimensionality reduction --------------------------

### 4.4.1. UMAP by condition (plot not included in the paper)
seu <- lapply(seu, function(seuObj){
  condition <- as.character(seuObj@meta.data$hash.ID)
  condition[condition == 'hto-mm'] <- 'mismatch control'
  condition[condition == 'hto-kd'] <- 'RUNX1::RUNX1T1 knockdown'
  seuObj@meta.data <- cbind(seuObj@meta.data, condition)
  return(seuObj)
})

p <- map2(seu, names(seu),
          ~ DimPlot(.x, 
                    reduction = "umap",
                    cols = cbPalette,
                    group.by = "condition") +
            theme_void() +
            labs(title = .y))
pdf(paste0(wd, '210_umap_cond.pdf'), height = 3, width = 8)
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'bottom')
dev.off()

### 4.4.2. UMAP by SoupOrCell assignment (Supplementary Figure 4B)
pdf(paste0(wd, "220_umap_soup.pdf"), height = 3, width = 8)

p <- map2(list(seu[[1]], seu[[2]]), names(seu)[1:2],
          ~DimPlot(.x,
                   group.by = "soup_assignment",
                   cols = c(cbPalette[3], cbPalette[5])) +
            theme_void() +
            labs(title = .y,
                 col = 'donor (SNP inference)'))

p3 <- DimPlot(seu[[3]],
              group.by = "soup_assignment",
              cols = c(cbPalette[5], cbPalette[3])) +
  theme_void() +
  labs(title = names(seu)[3],
       col = 'donor (SNP inference)')

l <- get_legend(p3)  

ggarrange(p[[1]], p[[2]], p3,
          ncol = 3,
          common.legend = TRUE,
          legend.grob = l,
          legend = 'bottom')
dev.off()

### 4.4.3. UMAP by AML and MSC markers (Supplementary Figure 4C)
p <- map2(seu, names(seu),
          ~ FeaturePlot(.x, 
                        features = c("LYZ", "NT5E", "THY1"),
                        ncol = 4) +
            plot_annotation(title = .y) &
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank()) 
)
pdf(paste0(wd, "225_umap_AML_MSCmarkers.pdf"), height = 7, width = 9)
ggarrange(plotlist = p, nrow = 3)
dev.off() #ITGAM = CD11b, FNT5E = CD73, THY1 = CD90

# ---

### 4.4.4. UMAP by cell cycle stage (Supplementary Figure 4E)
p <- map2(seu, names(seu),
          ~ DimPlot(.x,
                    reduction = "umap",
                    cols = cbPalette,
                    group.by = "Phase") +
            theme_void() +
            labs(title = .y,
                 col = 'Inferred cell cycle stage'))
pdf(paste0(wd, '230_umap_cellCycle.pdf'), height = 3, width = 8)
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'bottom')
dev.off()

# 5. Filter out MSCs ===========================================================

## 5.1. Throw away MSCs --------------------------------------------------------
seu[['patientA']] <- subset(seu[['patientA']], soup_assignment == '0')
seu[['patientB']] <- subset(seu[['patientB']], soup_assignment == '0')
seu[['patientC']] <- subset(seu[['patientC']], soup_assignment == '1')

## 5.2. Vizualize and make sure the correct cells were thrown away -------------

### 5.2.1. UMAP by SoupOrCell assignment (plot not included in the paper)
pdf(paste0(wd, "250_noMSC_umap_soup.pdf"), height = 3, width = 8)
p <- map2(list(seu[[1]], seu[[2]]), names(seu)[1:2],
          ~DimPlot(.x,
                   group.by = "soup_assignment",
                   cols = c(cbPalette[3], cbPalette[5])) +
            theme_void() +
            labs(title = .y,
                 col = 'donor (SNP inference)'))

p3 <- DimPlot(seu[[3]],
              group.by = "soup_assignment",
              cols = c(cbPalette[5], cbPalette[3])) +
  theme_void() +
  labs(title = names(seu)[3],
       col = 'donor (SNP inference)')

l <- get_legend(p3)  

ggarrange(p[[1]], p[[2]], p3,
          ncol = 3,
          common.legend = TRUE,
          legend.grob = l,
          legend = 'bottom')
dev.off()

# 6. Distinguish between leukemic and non-leukemic cells =======================

## 6.1. Cluster ----------------------------------------------------------------
seu[['patientA']] <- FindNeighbors(seu[['patientA']], dims = 1:35) 
seu[['patientB']] <- FindNeighbors(seu[['patientB']], dims = 1:29) 
seu[['patientC']] <- FindNeighbors(seu[['patientC']], dims = 1:19)

seu <- map(seu, ~ FindClusters(.x, resolution = 0.25, graph.name = 'SCT_snn'))

## 6.2. Vizualize clustering (plot not included in the paper) ------------------

p <- map2(seu, names(seu),
           ~ DimPlot(.x,
                     reduction = 'umap',
                     cols = cbPalette2,
                     group.by = 'SCT_snn_res.0.25',
                     label = TRUE) +
             plot_annotation(title = .y) +
             theme_void() +
             labs(title = ''))
pdf(paste0(wd, '260_noMSC_umap_cluster.pdf'), height = 3, width = 6)
ggarrange(plotlist = p, ncol = 3, legend = 'none')
dev.off()

## 6.3. Vizualize RUNX1T1 expression density (Figure 4A) -----------------------
p <- map(seu, 
          ~ plot_density(.x, 
                         reduction = "umap",
                         features = "RUNX1T1") +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank()) +
            labs(title = ''))
pdf(paste0(wd, '265_noMSC_umap_ETOdensity.pdf'), height = 3, width = 8)
ggarrange(plotlist = p, ncol = 3, common.legend = TRUE, legend = 'right')
dev.off()

## 6.4. Find cluster markers ---------------------------------------------------
myFindAllMarkers <- function(seuObj, ident) {
  Idents(seuObj) <- ident
  markers <- FindAllMarkers(seuObj, 
                            assay = "RNA",
                            slot = "data",
                            only.pos = FALSE, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25) %>%
    dplyr::group_by(cluster) 
  return(markers)
}

markers <- map(seu, ~ myFindAllMarkers(.x, 'SCT_snn_res.0.25'))
saveRDS(markers, paste0(wd, '270_SCTsnnRes0.25_markers.rds'))

##### get top 10 significant markers for each cluster
markers_filt <- map(markers, ~ dplyr::filter(.x, p_val_adj < 0.05) %>% 
                      arrange(desc(avg_log2FC)) %>%
                      dplyr::slice_max(avg_log2FC, n = 10))

## 6.5. Vizualize cluster markers (plot not included in the paper) -------------

pdf(paste0(wd, '280_heatmap_clust_markers.pdf'), height = 8, width = 5)
pmap(list(seu, names(seu), markers_filt),
     ~ DoHeatmap(..1, 
                 group.by = 'SCT_snn_res.0.25',
                 features = ..3$gene,
                 assay = "SCT",
                 slot = "scale.data",
                 group.bar = TRUE,
                 group.colors = cbPalette2) + 
       NoLegend() +
       scale_fill_viridis() +
       plot_annotation(title = ..2))
dev.off()


##### Patient A: clusters 5,6,7,9 are putative healthy cells because they have a low ETO density.
##### Patient B: clusters 5,7 are putative healthy cells because they have a low ETO density.
##### Patient C: all clusters are leukemic, which is to be expected, because this is PDX.

## 6.6. Make Seurat objects for healthy and leukemic cells, and also save the Seurat objects with both healthy and leukemic cells ---------------------

seu_healthy <- seu[c(1,2)]
seu_healthy[['patientA']] <- subset(seu_healthy[['patientA']], SCT_snn_res.0.25 == 5 | SCT_snn_res.0.25 == 6 | SCT_snn_res.0.25 == 7 | SCT_snn_res.0.25 == 9)
seu_healthy[['patientB']] <- subset(seu_healthy[['patientB']], SCT_snn_res.0.25 == 5 | SCT_snn_res.0.25 == 7)
saveRDS(seu_healthy, paste0(wd, '290_seu_healthy.rds'))

seu_leukemic <- seu
seu_leukemic [['patientA']] <- subset(seu_leukemic[['patientA']], SCT_snn_res.0.25 != 5 & SCT_snn_res.0.25 != 6 & SCT_snn_res.0.25 != 7 & SCT_snn_res.0.25 != 9)
seu_leukemic[['patientB']] <- subset(seu_leukemic[['patientB']], SCT_snn_res.0.25 != 5 & SCT_snn_res.0.25 != 7)
saveRDS(seu_leukemic, paste0(wd, '300_seu_leukemic.rds'))

saveRDS(seu, paste0(wd, '305_seu_complete.rds'))

# 7. Session info ==============================================================
sink(paste0(wd, '999_sessionInfo_1.txt'))
sessionInfo()
sink()

# 8. Clean the environment =====================================================
rm(list = ls())