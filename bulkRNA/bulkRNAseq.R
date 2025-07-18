library(tidyverse)
library(readxl)

##### differential expression analysis
library(DESeq2)

##### GSEA
library(msigdbr)
library(clusterProfiler)

##### vizualisations
library(ComplexHeatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(circlize)
library(viridis)
library(enrichplot)
library(patchwork)
library(aplot)
library(enrichplot)
library(ggpubr)

##### Replace with your working directory
setwd('C:/polina/analysis/pderevianko/runx1-runx1t1-kd-aml/bulkRNA/')

##### output directory (inside the working directory)
out_dir <- 'out/'

set.seed(2023)

kdmm_palette <- c('KD' = '#f41626', 'MM'= '#2538a5')

##### color blind-friendly color palette
cbPalette2 <- colorBlindness::paletteMartin
names(cbPalette2) <- NULL
cbPalette2 <- c(cbPalette2[3:9], cbPalette2[11:15])

##### Custom function to make a NES barplot after GSEA
plotNES <- function(GSEAres) {
  require(stringi)
  require(ggpubr)
  tib_toPlot <- GSEAres@result
  tib_toPlot$Significant <- ifelse(tib_toPlot$p.adjust < 0.05, 'yes', 'no')
  tib_top <- tib_toPlot %>% dplyr::filter(NES > 0) %>% arrange(NES) 
  tib_bot <- tib_toPlot %>% dplyr::filter(NES < 0) %>% arrange(NES) 
  max_x <- c(tib_top$NES, abs(tib_bot$NES)) %>% max()
  width_dif <- max(stri_length(tib_top$ID)) - max(stri_length(tib_bot$ID))
  
  p_top <- ggplot(tib_top,
                  aes(y = reorder(ID, -NES), x = NES)) +
    geom_col(aes(fill = -log10(qvalue))) +
    scale_fill_gradient(low = 'blue', high ='red') +
    geom_col(color = ifelse(tib_top$Significant == 'yes', 'black', NA), fill = NA, linewidth = 1) +
    scale_y_discrete(position = "right") +
    xlim(c(0, max_x)) +
    labs(x = NULL,
         y = NULL) +
    theme(plot.title = element_text(hjust = 0)) +
    theme_minimal() 
  
  p_bot <- ggplot(tib_bot,
                  aes(y = reorder(ID, -NES), x = NES)) +
    geom_col(aes(fill = -log10(qvalue))) +
    scale_fill_gradient(low = 'blue', high ='red') +
    geom_col(color = ifelse(tib_bot$Significant == 'yes', 'black', NA), fill = NA, linewidth = 1) +
    scale_y_discrete(position = "left") +
    xlim(c(-max_x, 0)) +
    labs(x = NULL,
         y = NULL) +
    theme(plot.title = element_text(hjust = 0)) +
    theme_minimal() 
  
  if (width_dif > 0) {
    p_top <- p_top + theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1),
                           plot.margin = unit(c(0, 0, 0, 6), 'pt'))
    p_bot <- p_bot + theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1),
                           plot.margin = unit(c(0, 6, 0, width_dif*6), 'pt'))
  }
  if (width_dif < 0) {
    p_top <- p_top + theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1),
                           plot.margin = unit(c(0, -width_dif*6, 0, 6), 'pt'))
    p_bot <- p_bot + theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1),
                           plot.margin = unit(c(0, 6, 0, 0), 'pt'))
  } 
  if (width_dif == 0) {
    p_top <- p_top + theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1),
                           plot.margin = unit(c(0, 0, 0, 6), 'pt'))
    p_bot <- p_bot + theme(axis.text.y = element_text(angle = 0, vjust = 1, hjust = 1),
                           plot.margin = unit(c(0, 6, 0, 0), 'pt'))
  } 
  
  p <- ggarrange(p_bot, p_top,
                 ncol = 2, widths = c(1,1), common.legend = TRUE)
  
  p <- annotate_figure(p, left = text_grob('Pathway', rot = 90), bottom = text_grob('Normalized enrichment score'))
  
  return(p)
}

##### Custom function to find out which genes in which GSEA gene set have which log2FC and padj
getPathwayGenes <- function(difexprTibble, geneset) {
  difexprTibbleFilt <- dplyr::filter(difexprTibble, gene %in% geneset)
  return(difexprTibbleFilt)
}

# 1. Load in data ==============================================================

## 1.1. Function to read input count files -------------------------------------
read_counts <- function(file_path) {
  file_content <- readLines(file_path)
  ##### extract sample name from the first line of metadata
  samplename <- strsplit(file_content[1], ' ')[[1]][3]
  ##### Use the fourth line (now first line of data_content) as column names
  col_names <- strsplit(file_content[4], "\t")[[1]]
  col_names[1] <- sub('^# ', '', col_names[1])
  ##### Remove the first three lines of metadata
  data_content <- file_content[-c(1:4)]
  ##### Read the data into a dataframe
  df <- read.table(text = data_content, 
                     header = FALSE, 
                     sep = "\t", 
                     col.names = col_names, 
                     quote = "", 
                     comment.char = "")
  return(list('samplename' = samplename,
              'data' = df))
}

## 1.2. Read the count files to dataframes -------------------------------------
in_dir <- 'in'

##### These count files are available at https://doi.org/10.5281/zenodo.14578307, folder "bulkRNA"
file_paths <- list.files(in_dir, pattern = "\\.txt$", full.names = TRUE)

in_data <- lapply(file_paths, read_counts) 
in_data <- setNames(lapply(in_data, `[[`, 2), sapply(in_data, `[[`, 1))

# 2. Prepare the data for DESeq2 ===============================================

## 2.1. Extract raw counts. Sum them up for different isoforms of same genes ----
counts_list <- lapply(names(in_data), function(samplename){
  df_filt <- dplyr::select(in_data[[samplename]], Counts, GeneName)
  df_aggr <- aggregate(Counts ~ GeneName,
                       data = df_filt,
                       FUN = sum)
  colnames(df_aggr)[colnames(df_aggr) == 'Counts'] <- samplename
  return(df_aggr)
})
names(counts_list) <- names(in_data)

## 2.2. Make count matrix ------------------------------------------------------
counts <- dplyr::select(counts_list[[1]], GeneName)
for (samplename in names(in_data)){
  counts <- left_join(counts, counts_list[[samplename]], by = 'GeneName')
}

rownames(counts) <- counts$GeneName
counts$GeneName <- NULL

## 2.3. Make metadata dataframe ------------------------------------------------
coldata <- data.frame('condition' = as.factor(c(rep('KD',3), rep('MM',3))))
rownames(coldata) <- colnames(counts)

# 3. DESeq2 ====================================================================

## 3.1. Run --------------------------------------------------------------------
deseq_data <- DESeqDataSetFromMatrix(countData = counts,
                                     colData = coldata,
                                     design = ~ condition)
deseq_data$condition <- factor(deseq_data$condition, levels = c('MM', 'KD'))
dds <- DESeq(deseq_data)

## 3.2. Save results (not provided in the paper materials) ---------------------
saveRDS(dds, paste0(out_dir, '015_deseq2_out.rds'))
##### dds <- readRDS(paste0(out_dir, '015_deseq2_out.rds'))

# 4. Transformations and filtering of DESeq2 output ============================

## 4.1. Make a data frame with DESeq2 output (Suppl. Table 6) ------------------
res <- results(dds)

res_df <- as.data.frame(res) %>% 
  merge(counts, by = 'row.names', all = TRUE)
rownames(res_df) <- res_df$Row.names
res_df$Row.names <- NULL

write.csv(res_df, paste0(out_dir, '020_deseq2_out.csv'))

## 4.2. Transform --------------------------------------------------------------
rld <- rlog(dds, blind = FALSE)

## 4.3. Calculate z-scores -----------------------------------------------------
z <- t(apply(assay(rld), 1, scale)) 
colnames(z) <- colnames(rld)

## 4.4. Remove NAs from the results --------------------------------------------
res_ex <- drop_na(res_df, log2FoldChange, padj)

## 4.5. Filter by significance and log2FC --------------------------------------
res_sign <- dplyr::filter(res_ex, padj < 0.05 & baseMean > 50)

res_sel_top <- slice_max(res_sign, log2FoldChange, n = 50) %>% arrange(desc(log2FoldChange))
res_sel_bot <- slice_min(res_sign, log2FoldChange, n = 50) %>% arrange(desc(log2FoldChange))
res_sel <- bind_rows(res_sel_top, res_sel_bot)
res_sel$gene <- rownames(res_sel)

z_sel <- z[rownames(res_sel),] 

l2_val <- as.matrix(res_sel$log2FoldChange)
colnames(l2_val) <- 'log2FC'

means <- as.matrix(res_sel$baseMean)
colnames(means) <- 'AveExpr'

# 5. Visualize DESeq2 results ==================================================

## 5.1. PCA plot for differences between replicates (Suppl. Figure 3A) ---------
pdf(paste0(out_dir, '030_rld_PCA.pdf'), height = 5, width = 5)
plotPCA(rld, intgroup = 'condition') +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = kdmm_palette, 
                     labels = c('mismatch control', 'RUNX1::RUNX1T1 knockdown'), 
                     name = 'Condition') +
  ylim(-20, 20) +
  theme(legend.position = 'bottom')
dev.off()

## 5.2. Heatmap of z-scores for top and bottom 50 DE-genes (Figure 3B) ---------

### 5.2.1. Color palettes
col_z <- colorRamp2(seq(min(z_sel), max(z_sel), length.out = nrow(z_sel)), 
                        viridis(nrow(z_sel)))

col_log2fc <- colorRamp2(seq(min(l2_val), max(l2_val), length.out = nrow(l2_val)),
                         viridis(nrow(l2_val)))

col_AveExpr <- colorRamp2(seq(-quantile(means, probs = 0.9), quantile(means, probs = 0.9), length.out = nrow(means)),
                          viridis(nrow(means)))

### 5.2.2. Heatmap
pdf(paste0(out_dir, '040_rld-z_top50bottom50DEgenes.pdf'), width = 6, height = 18)
h1 <- Heatmap(z_sel, 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_labels = colnames(z_sel),
              name = 'Z-score',
              col = col_z,
              width = unit(1.5, 'inches'),
              top_annotation = HeatmapAnnotation(condition = coldata$condition,
                                                 col = list(condition = kdmm_palette)))
h2 <- Heatmap(l2_val,
              row_labels = rownames(z_sel),
              cluster_rows = FALSE,
              name = 'log2FC',
              col = col_log2fc,
              width = unit(0.7, 'inches'),
              cell_fun = function(j, i, x, y, w, h, col){
                grid.text(round(l2_val[i, j], 2), x, y)
              })
h3 <- Heatmap(means,
              row_labels = rownames(z_sel),
              cluster_rows = FALSE,
              name = 'AveExpr',
              col = col_AveExpr,
              width = unit(0.7, 'inches'),
              row_names_gp = gpar(fontface = "bold"),
              cell_fun = function(j, i, x, y, w, h, col){
                grid.text(round(means[i, j], 2), x, y)
              })
h1+h2+h3
dev.off()

## 5.3. Volcano plots ----------------------------------------------------------

### 5.3.1. Volcano plot with all default genes labelled (Suppl. Figure 3D)
pdf(paste0(out_dir, '050_volcano_allGenesLabeled.pdf'),
    width = 8,
    height = 8)
EnhancedVolcano(res_ex,
                lab = rownames(res_ex),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-3,
                FCcutoff = 2,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                boxedLabels = FALSE,
                col = c('#444444','#444444', cbPalette2[7], cbPalette2[2]),
                colAlpha = 1,
                pointSize = 0.5,
                labSize = 2,
                title = NULL,
                subtitle = NULL
)
dev.off()

### 5.3.2. Volcano plot with target genes of RUNX1::RUNX1T1 and differentiation markers (Figure 3C)
genes_of_interest <- c('CD34',
                       'LINC01257',
                       'ANGPT1',
                       'CEBPA',
                       'CEBPE',
                       'CLEC12A',
                       'C11orf21',
                       'EPX',
                       'RNASE2',
                       'RNASE3',
                       'RUNX1T1',
                       'OGG1',
                       'CEACAM6',
                       'CAS6',
                       'NKG7',
                       'IGFBP7',
                       'PRKCD',
                       'LYZ',
                       'ELANE',
                       'CD24',
                       'PLAC8',
                       'CTSG',
                       'PRG2',
                       'PRG3',
                       'CLC',
                       'ABCA13',
                       'CYBB',
                       'DNAH8',
                       'L1CAM',
                       'RHOXF1P1',
                       'TGM5',
                       'MS4A3',
                       'ADAM12',
                       'CST7',
                       'NPTX1',
                       'IL32',
                       'ACSM3',
                       'MMP1',
                       'CNTNAP5',
                       'KRT15',
                       'BPI',
                       'CEACAM8',
                       'TSPAN32',
                       'PCBP3',
                       'FCRL1',
                       'FCRL2',
                       'MYBPH',
                       'CA9',
                       'FUT7',
                       'MPO',
                       'VSIR')

pdf(paste0(out_dir, '060_volcano_REtargetGenes_diffMarkers.pdf'),
    width = 8,
    height = 8)
EnhancedVolcano(res_ex,
                lab = rownames(res_ex),
                selectLab = genes_of_interest,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 10e-3,
                FCcutoff = 1,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                boxedLabels = TRUE,
                col = c('#444444','#444444', cbPalette2[7], cbPalette2[2]),
                colAlpha = 1,
                pointSize = 1,
                labSize = 4,
                title = NULL,
                subtitle = NULL
)
dev.off()

# 6. GSEA ======================================================================

## 6.1. Get the list of gene sets we would like to use for analysis ------------

mygenesets_names <- c(
  "TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_HSC_DN",
  "TONKS_TARGETS_OF_RUNX1_RUNX1T1_FUSION_HSC_UP",
  "HAY_BONE_MARROW_CD34_POS_CLP",
  "HAY_BONE_MARROW_CD34_POS_EO_B_MAST",
  "HAY_BONE_MARROW_CD34_POS_ERP",
  "HAY_BONE_MARROW_CD34_POS_ERP_EARLY",
  "HAY_BONE_MARROW_CD34_POS_GRAN",
  "HAY_BONE_MARROW_CD34_POS_HSC",
  "HAY_BONE_MARROW_CD34_POS_LMPP",
  "HAY_BONE_MARROW_CD34_POS_LYMPHOID_UNK",
  "HAY_BONE_MARROW_CD34_POS_MEP",
  "HAY_BONE_MARROW_CD34_POS_MKP",
  "HAY_BONE_MARROW_CD34_POS_MULTILIN",
  "HAY_BONE_MARROW_CD34_POS_PRE_B",
  "HAY_BONE_MARROW_CD34_POS_PRE_PC",
  "HAY_BONE_MARROW_CD8_T_CELL",
  "HAY_BONE_MARROW_DENDRITIC_CELL",
  "HAY_BONE_MARROW_EARLY_ERYTHROBLAST",
  "HAY_BONE_MARROW_EOSINOPHIL",
  "HAY_BONE_MARROW_ERYTHROBLAST",
  "HAY_BONE_MARROW_FOLLICULAR_B_CELL",
  "HAY_BONE_MARROW_IMMATURE_NEUTROPHIL",
  "HAY_BONE_MARROW_MONOCYTE",
  "HAY_BONE_MARROW_NAIVE_T_CELL",
  "HAY_BONE_MARROW_NEUTROPHIL",
  "HAY_BONE_MARROW_NK_CELLS",
  "HAY_BONE_MARROW_PLASMA_CELL",
  "HAY_BONE_MARROW_PLATELET",
  "HAY_BONE_MARROW_PRE_DENDRITIC",
  "HAY_BONE_MARROW_PRO_B",
  "HAY_BONE_MARROW_STROMAL"
)

mygenesets <- msigdbr(species = 'Homo sapiens') %>%
  dplyr::filter(gs_name %in% mygenesets_names)

mygenesets_list <- split(mygenesets, 
                         x = mygenesets$gene_symbol, 
                         f = mygenesets$gs_name)

## 6.2. Make gene ranks --------------------------------------------------------

##### I use the stat column of the DESeq2 result as the ranking metric.

ranks <- res_df$stat
names(ranks) <- rownames(res_df)
ranks <- ranks[!is.na(ranks)] %>% 
  sort(decreasing = TRUE) ##### sorting is required for ClusterProfiler

## 6.3. Run GSEA ---------------------------------------------------------------
gsea_res <- GSEA(
  geneList = ranks,
  maxGSSize = 2000,
  pvalueCutoff = 1,
  eps = 0,
  pAdjustMethod = 'BH',
  seed = TRUE,
  TERM2GENE = dplyr::select(mygenesets, gs_name, gene_symbol)
)

saveRDS(gsea_res, paste0(out_dir, '071_gsea_res.rds'))
##### gsea_res <- readRDS(paste0(out_dir, '071_gsea_res.rds'))

gsea_sign <- dplyr::filter(gsea_res, p.adjust < 0.05)

## 6.4. Plot GSEA-style plots for all the analyzed pathways (Suppl. Figure 3E for 3 of these pathways) ----
pdf(paste0(out_dir, '080_gseaPlots.pdf'), width = 5, height = 5)
plotlist <- map(mygenesets_names[mygenesets_names %in% gsea_sign$ID],
    function(geneset_name){
      p <- enrichplot::gseaplot(
        gsea_sign,
        geneSetID = geneset_name,
        title = geneset_name) %>%
        as.patchwork()
      p <- p + plot_annotation(subtitle = paste0(
        'NES = ', as.character(signif(gsea_res@result$NES[gsea_res@result$ID == geneset_name], digits = 3)),
        ', padj = ', as.character(signif(gsea_res@result$p.adjust[gsea_res@result$ID == geneset_name], digits = 3)),
        ', FDR = ', as.character(signif(gsea_res@result$qvalue[gsea_res@result$ID == geneset_name], digits = 3))
      ))
      return(p)
    })
print(plotlist)
dev.off()

## 6.5. Visualize GSEA results for selected pathways ---------------------------

### 6.5.1. Determine if this set of gene sets has overlapping genes
hay_names <- mygenesets_names[grep('^HAY_', mygenesets_names)]
hay <- mygenesets_list[hay_names]

hay_allGenes <- unlist(hay)
hay_geneCounts <- table(hay_allGenes)
names(hay_geneCounts[hay_geneCounts > 1])
##### The Hay set of gene sets only has one repeating gene (CYB561D2) so I would like to make a dot/barplot with it.

gsea_hay <- dplyr::filter(gsea_res, grepl('^HAY_', ID))

#### 6.5.2. Barplot (Figure 3D)
pdf(paste0(out_dir, '100_gsea_HAY_BONE_MARROW_barplot.pdf'), width = 10, height = 6)
plotNES(gsea_hay)
dev.off()

# 7. Compare with the RNAseq results of Kasumi-1 and SKNO-1 ====================

## 7.1. Read in the RNAseq results of Kasumi-1 and SKNO-1 ----------------------

##### These input files are from Issa et al 2023 and are available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217113
cellines_df <- read_excel('in/LNP_K1_S1_RNAseq.xlsx', sheet = 'All')

cellines_samplenames <- c('K1_mm_rep1', 'K1_kd_rep1', 
                          'K1_mm_rep2', 'K1_kd_rep2',
                          'K1_mm_rep3', 'K1_kd_rep3',
                          'S1_mm_rep1', 'S1_kd_rep1', 
                          'S1_mm_rep2', 'S1_kd_rep2',
                          'S1_mm_rep3', 'S1_kd_rep3')

kasumi_samplenames <- cellines_samplenames[grepl('^K1_', cellines_samplenames)]
skno_samplenames <- cellines_samplenames[grepl('^S1_', cellines_samplenames)]

colnames(cellines_df)[9:20] <- cellines_samplenames

## 7.2. Extract the counts for each sample separately --------------------------
cellines_counts <- map(cellines_samplenames, 
                       ~ dplyr::select(cellines_df, GeneName, !!sym(.x)))
names(cellines_counts) <- cellines_samplenames

## 7.3 Sum up the counts for different isoforms of same genes ------------------
cellines_counts <- lapply(names(cellines_counts), function(samplename){
    df <- cellines_counts[[samplename]]
    colnames(df)[2] <- 'Counts'
    df_aggr <- aggregate(Counts ~ GeneName,
                         data = df,
                         FUN = sum)
    colnames(df_aggr)[colnames(df_aggr) == 'Counts'] <- samplename
    return(df_aggr)
  })
names(cellines_counts) <- cellines_samplenames

## 7.4. Make count matrices ----------------------------------------------------
kasumi_counts <- cellines_counts[kasumi_samplenames]
skno_counts <- cellines_counts[skno_samplenames]
cellines_counts <- list('Kasumi-1' = kasumi_counts, 'SKNO-1' = skno_counts)

cellines_mtx <- lapply(cellines_counts, function(df_list){
  mtx <- dplyr::select(df_list[[1]], GeneName)
  for(samplename in names(df_list)){
    mtx <- left_join(mtx, df_list[[samplename]], by = 'GeneName')
  }
  rownames(mtx) <- mtx$GeneName
  mtx$GeneName <- NULL
  return(mtx)
})

## 7.5. Make metadata dataframes -----------------------------------------------
coldata <- data.frame('condition' = as.factor(rep(c('MM', 'KD'),3)))
coldata <- list('Kasumi-1' = coldata, 'SKNO-1' = coldata)
coldata <- lapply(names(coldata), function(celline_name){
  rownames(coldata[[celline_name]]) <- colnames(cellines_mtx[[celline_name]])
  return(coldata[[celline_name]])
})
names(coldata) <- names(cellines_mtx)

## 7.6. DESeq2 -----------------------------------------------------------------
deseq_cellines <- lapply(names(cellines_mtx), function(celline_name){
  deseq_dat <- DESeqDataSetFromMatrix(countData = cellines_mtx[[celline_name]],
                                      colData = coldata[[celline_name]],
                                      design = ~ condition)
  deseq_dat$condition <- factor(deseq_dat$condition, levels = c('MM', 'KD'))
  dds_celline <- DESeq(deseq_dat)
  return(dds_celline)
})
names(deseq_cellines) <- names(cellines_mtx)

## 7.7. Save results (not provided in the paper materials) ---------------------
saveRDS(deseq_cellines, paste0(out_dir, '110_cellines_deseq2_out.rds'))
##### deseq_cellines <- readRDS(paste0(out_dir, '110_cellines_deseq2_out.rds'))

## 7.8. Transformations and filtering of DESeq2 output -------------------------
res_cellines <- map(deseq_cellines, results)

res_cellines_df <- lapply(names(res_cellines), function(celline_name){
  df <- as.data.frame(res_cellines[[celline_name]]) %>%
    merge(cellines_mtx[[celline_name]], by = 'row.names', all = TRUE)
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  return(df)
})
names(res_cellines_df) <- names(res_cellines)

res_ex_cellines <- map(res_cellines_df, ~ drop_na(.x, log2FoldChange, padj))

## 7.9. Filter the DESeq2 output for the cell lines for significance and log2FC ----
res_cellines_filt <- lapply(res_ex_cellines, function(df){
  df_new <- df %>% dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 1)
  return(df_new)
})

## 7.10. Split to down- and upregulated genes ----------------------------------
res_cellines_split <- lapply(res_cellines_filt, function(df){
  sublist <- list(
    'up' = dplyr::filter(df, log2FoldChange > 0),
    'dn' = dplyr::filter(df, log2FoldChange < 0)
  )
})

## 7.11. Extract the names of the genes that are substantially & significantly differentially expressed in each of the cell lines ----
cellines_genesets <- lapply(res_cellines_split, function(sublist){
  sublist_new <- lapply(sublist, function(df){
    vec <- rownames(df)
    return(vec)
  })
  return(sublist_new)
})

cellines_genesets <- lapply(names(cellines_genesets), function(celline_name){
  sublist <- lapply(names(cellines_genesets[[celline_name]]), function(up_or_dn){
    gene_symbol <- cellines_genesets[[celline_name]][[up_or_dn]]
    gs_name <- paste0(rep(celline_name, length(gene_symbol)), '_', up_or_dn)
    df <- data.frame(gs_name, gene_symbol)
    return(df)
  })
  return(sublist)
})

cellines_genesets <- lapply(cellines_genesets, function(sublist){
  df <- rbind(sublist[[1]], sublist[[2]])
  return(df)
}) 
names(cellines_genesets) <- names(res_cellines_filt)

## 7.12. Save results (not provided in the paper materials) --------------------
saveRDS(cellines_genesets, paste0(out_dir, '140_cellines_genesets.rds'))
##### cellines_genesets <- readRDS(paste0(out_dir, '140_cellines_genesets.rds'))

## 7.13. Run GSEA of the PDX against the genesets from cellines ----------------
gsea_res_cellines <- map(cellines_genesets, ~GSEA(
  geneList = ranks,
  maxGSSize = 2000,
  pvalueCutoff = 1,
  eps = 0,
  pAdjustMethod = 'BH',
  seed = TRUE,
  TERM2GENE = .x
))

## 7.14. Save GSEA results (not provided in paper materials) -------------------
saveRDS(gsea_res_cellines, paste0(out_dir, '151_gsea_res_againstCellines.rds'))
##### gsea_res_cellines <- readRDS(paste0(out_dir, '151_gsea_res_againstCellines.rds'))

## 7.15. Make GSEA-style plots (Figure 3A, Suppl. Figures 3B–C) ----------------
pdf(paste0(out_dir, '160_gseaPlots_againstCellines.pdf'), width = 11, height = 11)
plotlist1 <- map(c('Kasumi-1_up', 'Kasumi-1_dn'),
                function(celline_name){
                  p <- enrichplot::gseaplot(
                    gsea_res_cellines[['Kasumi-1']],
                    geneSetID = celline_name,
                    title = celline_name,
                    color.line = '#21908c') %>%
                    as.patchwork()
                  p <- p + plot_annotation(subtitle = paste0(
                    'NES = ', as.character(signif(gsea_res_cellines[['Kasumi-1']]@result$NES[gsea_res_cellines[['Kasumi-1']]@result$ID == celline_name], digits = 3)),
                    ', padj = ', as.character(signif(gsea_res_cellines[['Kasumi-1']]@result$p.adjust[gsea_res_cellines[['Kasumi-1']]@result$ID == celline_name], digits = 3)),
                    ', FDR = ', as.character(signif(gsea_res_cellines[['Kasumi-1']]@result$qvalue[gsea_res_cellines[['Kasumi-1']]@result$ID == celline_name], digits = 3))
                  ))
                  return(p)
                })
plotlist2 <- map(c('SKNO-1_up', 'SKNO-1_dn'),
                 function(celline_name){
                   p <- enrichplot::gseaplot(
                     gsea_res_cellines[['SKNO-1']],
                     geneSetID = celline_name,
                     title = celline_name,
                     color.line = '#21908c') %>%
                     as.patchwork()
                   p <- p + plot_annotation(subtitle = paste0(
                     'NES = ', as.character(signif(gsea_res_cellines[['SKNO-1']]@result$NES[gsea_res_cellines[['SKNO-1']]@result$ID == celline_name], digits = 3)),
                     ', padj = ', as.character(signif(gsea_res_cellines[['SKNO-1']]@result$p.adjust[gsea_res_cellines[['SKNO-1']]@result$ID == celline_name], digits = 3)),
                     ', FDR = ', as.character(signif(gsea_res_cellines[['SKNO-1']]@result$qvalue[gsea_res_cellines[['SKNO-1']]@result$ID == celline_name], digits = 3))
                   ))
                   return(p)
                 })
print(plotlist1)
print(plotlist2)
dev.off()

plotlist <- append(plotlist1, plotlist2)

svg(paste0(out_dir, '160_gseaPlots_againstCellines.svg'))
ggarrange(plotlist = plotlist, ncol = 2, nrow = 2)
dev.off()

# 99. Session info =============================================================
sink(paste0(out_dir, '999_sessionInfo.txt'))
sessionInfo()
sink()