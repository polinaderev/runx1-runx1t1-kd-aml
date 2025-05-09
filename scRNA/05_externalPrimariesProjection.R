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
wd <- 'repos/runx1eto-kd-aml/scRNA/out/'

##### reproducibility
set.seed(42)

# 1. Load in the Zeng et al reference ==========================================
##### Reference paper: https://www.biorxiv.org/content/10.1101/2023.12.26.573390v1
##### Reference Seurat object provided by Andy Zeng
ref <- readRDS('references/Zeng_BoneMarrowMap_Annotated_Dataset.rds')
Idents(ref) <- 'Tissue'

ref_data <- Embeddings(subset(ref, downsample = 50000)[["umap"]]) %>%
  as.data.frame() 

# 2. RUNX1::RUNX1T1 AMLs from Lambo et al 2023 =================================

## 2.1. Load in the data -------------------------------------------------------
##### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235063
ref_path <- 'references/Lambo2023/'

sample_names_dx <- paste0(
  c('GSM7494257_AML16_',
    'GSM7494266_AML15_',
    'GSM7494314_AML14_',
    'GSM7494329_AML13_',
    'GSM7494326_AML12_'),
  'DX_processed_'
)

sample_names_rel <- paste0(
  c('GSM7494258_AML16_',
    'GSM7494267_AML15_',
    ##### patient 14 was censored so we don't have their relapse dataset
    'GSM7494330_AML13_',
    'GSM7494327_AML12_'),
  'REL_processed_'
)

sample_names <- list(
  'dx' = sample_names_dx,
  'rel' = sample_names_rel
)

seu <- lapply(sample_names, function(sublist){
  sublist_new <- lapply(sublist, function(sample_name){
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
  return(sublist_new)
})

names(seu[['dx']]) <- c('AML16_DX', 'AML15_DX', 'AML14_DX', 'AML13_DX', 'AML12_DX')
names(seu[['rel']]) <- c('AML16_DX', 'AML15_DX', 'AML13_DX', 'AML12_DX')

## 2.2. Subset to malignant cells only -----------------------------------------
seu <- lapply(seu, function(sublist){
  sublist_new <- lapply(sublist, function(seuObj){
    seuObj_subset <- subset(seuObj, Malignant == 'Malignant')
    return(seuObj_subset)
  })
  return(sublist_new)
})

## 2.3. Project onto Andy Zeng's reference -------------------------------------
seu <- lapply(seu, function(sublist) {
  sublist_new <- lapply(sublist, function(seuObj){
  seuObj <- NormalizeData(seuObj)
  seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 2000)
  DefaultAssay(seuObj) <- 'RNA'
  return(seuObj)
  })
  return(sublist_new)
})

anchors <- lapply(seu, function(sublist){
  map(sublist, ~ FindTransferAnchors(reference = ref, 
                                     query = .x, 
                                     dims = 1:30, 
                                     npcs = 30))
  })

seu <- map2(anchors, seu, function(anchor_sublist, seu_sublist){
  res <- map2(anchor_sublist, seu_sublist,
       ~ MapQuery(anchorset = .x, 
                  reference = ref, 
                  query = .y, 
                  refdata = list(Zeng_celltype = 'CellType_Broad'),
                  reference.reduction = "pca", 
                  reduction.model = "umap"))
  return(res)
})

preds <- lapply(seu, function(sublist){
  sublist_new <- lapply(sublist, function(seuObj){
    df <- seuObj@meta.data %>%
      rownames_to_column('cell_label') %>% 
      select(cell_label, predicted.Zeng_celltype, predicted.Zeng_celltype.score)
    return(df)
  })
  return(sublist_new)
})

saveRDS(preds, paste0(wd, '701_zeng_lamboPreds.rds'))
# preds <- readRDS(paste0(wd, '701_zeng_lamboPreds.rds'))
#
# seu <- lapply(names(seu), function(diagnosis_or_relapse){
#   sublist_new <- lapply(names(seu[[diagnosis_or_relapse]]), function(patient){
#     seu[[diagnosis_or_relapse]][[patient]]@meta.data$cell_label <- rownames(seu[[diagnosis_or_relapse]][[patient]]@meta.data)
#     seu[[diagnosis_or_relapse]][[patient]]@meta.data <- seu[[diagnosis_or_relapse]][[patient]]@meta.data %>%
#       left_join(preds[[diagnosis_or_relapse]][[cell_label]] %>% rownames_to_column('cell_label'),
#                 by = 'cell_label') %>%
#       column_to_rownames(var = 'cell_label')
#     return(seu[[diagnosis_or_relapse]][[patient]])
#   })
#   names(sublist_new) <- names(seu[[diagnosis_or_relapse]])
#   return(sublist_new)
# })
# names(seu) <- names(preds)

pred_embeddings <- lapply(seu, function(sublist){
  sublist_new <- lapply(sublist, function(seuObj){
    pca <- Embeddings(seuObj[["ref.pca"]])
    umap <- Embeddings(seuObj[["ref.umap"]])
    out <- list('pca' = pca, 'umap' = umap)
    return(out)
  })
  return(sublist_new)
})

saveRDS(pred_embeddings, paste0(wd, '702_zeng_lamboPred_embeddings.rds'))
# pred_embeddings <- readRDS(paste0(wd, '702_zeng_lamboPred_embeddings.rds'))
# 
# seu <- map2(seu, pred_embeddings, function(seuObj, embList){
#   seuObj[['ref.pca']] <- CreateDimReducObject(as.matrix(embList[['pca']]))
#   seuObj[['ref.umap']] <- CreateDimReducObject(as.matrix(embList[['umap']]))
#   return(seuObj)
# })

## 2.3. Visualize --------------------------------------------------------------
query_data <- lapply(names(seu), function(diagnosis_or_relapse){
  sublist_new <- lapply(names(seu[[diagnosis_or_relapse]]), function(sample_name){
    seu[[diagnosis_or_relapse]][[sample_name]]@meta.data$cell_label <- rownames(seu[[diagnosis_or_relapse]][[sample_name]]@meta.data)
    df <- Embeddings(seu[[diagnosis_or_relapse]][[sample_name]][['ref.umap']]) %>%
      as.data.frame() %>%
      mutate(patient = sample_name) %>%
      rename(UMAP_1 = refUMAP_1,
             UMAP_2 = refUMAP_2) %>%
      tibble::rownames_to_column(var = 'cell_label') %>%
      left_join(seu[[diagnosis_or_relapse]][[sample_name]]@meta.data %>% 
                  select(cell_label, predicted.Zeng_celltype.score),
                by = 'cell_label')
    return(df)
  })
  names(sublist_new) <- names(seu[[diagnosis_or_relapse]])
  return(sublist_new)
})
names(query_data) <- names(seu)

### 2.3.1. UMAP colored by cell density
plotlist <- lapply(names(seu), function(diagnosis_or_relapse){
  plts <- lapply(names(seu[[diagnosis_or_relapse]]), function(sample_name){
    plt <- query_data[[diagnosis_or_relapse]][[sample_name]] %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(data = ref_data, color = '#E3E3E3', size = 0.05, alpha = 0.5) +
      ggpointdensity::geom_pointdensity(size = 0.2) +
      jcolors::scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 1.75) +
      geom_density_2d(alpha = 0.4, color = 'black', h = 1.5, linewidth = 0.3) +
      theme_void() + 
      labs(subtitle = sample_name) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 18), 
                     legend.position = 'none')
    return(plt)
  })
  names(plts) <- names(seu[[diagnosis_or_relapse]])
  return(plts)
})
names(plotlist) <- names(seu)

pdf(paste0(wd, '705_lambo_umapAndy_density.pdf'), height = 5.5, width = 14)

ggarrange(plotlist[['dx']][['AML12_DX']], plotlist[['rel']][['AML12_REL']], nrow = 1, ncol = 2) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML13_DX']], plotlist[['rel']][['AML13_REL']], nrow = 1, ncol = 2) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML14_DX']], ggplot() + theme_void(), nrow = 1, ncol = 2) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML15_DX']], plotlist[['rel']][['AML15_REL']], nrow = 1, ncol = 2) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML16_DX']], plotlist[['rel']][['AML16_REL']], nrow = 1, ncol = 2) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

dev.off()

### 2.3.3. UMAP colored by cell type prediction confidence

plotlist <- lapply(names(seu), function(diagnosis_or_relapse){
  plts <- lapply(names(seu[[diagnosis_or_relapse]]), function(sample_name){
    plt <- query_data[[diagnosis_or_relapse]][[sample_name]] %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(data = ref_data, color = '#E3E3E3', size = 0.05, alpha = 0.5) +
      geom_point(aes(colour = predicted.Zeng_celltype.score),
                 data = query_data[[diagnosis_or_relapse]][[sample_name]],
                 size = 0.05) +
      scale_color_viridis_c(begin = 0,
                            end = 1,
                            direction = 1,
                            na.value = 'grey',
                            name = 'Prediction confidence') +
      theme_void() + 
      labs(subtitle = sample_name) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 18), 
                     legend.position = 'bottom')
    return(plt)
  })
  names(plts) <- names(seu[[diagnosis_or_relapse]])
  return(plts)
})
names(plotlist) <- names(seu)

pdf(paste0(wd, '710_lambo_umapAndy_predConfidence.pdf'), height = 6, width = 14)

ggarrange(plotlist[['dx']][['AML12_DX']], plotlist[['rel']][['AML12_REL']], nrow = 1, ncol = 2, common.legend = FALSE) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML13_DX']], plotlist[['rel']][['AML13_REL']], nrow = 1, ncol = 2, common.legend = FALSE) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML14_DX']], ggplot() + theme_void(), nrow = 1, ncol = 2) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML15_DX']], plotlist[['rel']][['AML15_REL']], nrow = 1, ncol = 2, common.legend = FALSE) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

ggarrange(plotlist[['dx']][['AML16_DX']], plotlist[['rel']][['AML16_REL']], nrow = 1, ncol = 2, common.legend = FALSE) +
  plot_annotation(title = 'Lambo et al 2023, projected onto Zeng et al 2023')

dev.off()

# 98. Session info ==============================================================
sink(paste0(wd, '999_sessionInfo_3.txt'))
sessionInfo()
sink()