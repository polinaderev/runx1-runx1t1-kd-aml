library(flowCore)
library(Rtsne)
library(tidyverse)
library(readxl)
library(flownalysis)
library(FactoMineR)
library(CytoML)

##### vizualization
library(ggcyto)
library(ggpubr)

kdmm_palette <- c('RUNX1::RUNX1T1 knockdown' = '#f41626', 'mismatch control'= '#2538a5')
cbPalette <- c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#0072B2","#D55E00", "#CC79A7", "#999999")

in_path <- 'flow_ptABD/in/' ##### Folder with input files. Replace with your path
names(in_path) <- c('EX68', 'EX77')
wd <- 'flow_ptABD/out/' ##### output directory. Replace with your path

##### reproducibility
set.seed(42)

# 1. Load in the data ==========================================================

## 1.1. Metadata Excel file that contains info about samples and about the panel -----
##### Available at https://doi.org/10.5281/zenodo.14578307, folder "flow/flow_ptABD_gated"
excel_path <- paste0(in_path, 'metadata.xlsx')
metadata_sheets <- excel_sheets(excel_path)

metadata <- map(metadata_sheets, ~read_excel(excel_path, sheet = .x))
names(metadata) <- c('samples', 'panel')

## 1.2. FCS files --------------------------------------------------------------
##### Available at https://doi.org/10.5281/zenodo.14578307, folder "flow/flow_ptABD_gated"

### 1.2.1. Load 
dat_raw <- read.flowSet(files = metadata[['samples']]$file_name,
                        path = in_path,
                        descriptions = metadata[['samples']]$sample_id,
                        transformation = FALSE,
                        truncate_max_range = FALSE)

sampleNames(dat_raw) <- metadata[['samples']]$sample_id
pData(dat_raw)$name <- metadata[['samples']]$sample_id

### 1.2.2. Replace the column names (channels) with marker names
rename_vector <- setNames(metadata[['panel']]$marker,
                          metadata[['panel']]$channel)

columns_to_rename <- colnames(dat_raw) %in% names(rename_vector)

colnames(dat_raw)[columns_to_rename] <- rename_vector[colnames(dat_raw)[columns_to_rename]]

# 2. Prepare the data ==========================================================

## 2.1. Only protein marker columns will be used for the analysis --------------
marker_cols <- metadata[['panel']]$marker
dat_raw <- dat_raw[, marker_cols]

## 2.2. Remove cells with negative values --------------------------------------

### 2.2.1. Prepare gates for each marker of non-negative values
gates <- lapply(marker_cols, function(vec){
  vec_new <- map(vec,
                     function(marker){
                       mat <- matrix(c(0,Inf), ncol = 1, dimnames = list(c('min', 'max'), marker))
                       gate <- rectangleGate(filterId = paste0(sub('-.*', '', marker), '_clean'), .gate = mat)
                       return(gate)
                       })
  return(vec_new)
})

### 2.2.2. Gate and subset 
dat_filt <- lapply(names(in_path), function(experiment_name){
  dat_gated <- filter(dat_raw[[experiment_name]], gates[[experiment_name]][[1]])
  dat_new <- Subset(dat_raw[[experiment_name]], dat_gated)
  
  for (i in 2:length(gates[[experiment_name]])){
    dat_gated <- filter(dat_new, gates[[experiment_name]][[i]])
    dat_new <- Subset(dat_new, dat_gated)
  }
  
  return(dat_new)
})
names(dat_filt) <- names(in_path)

## 2.3. Arcsinh transformation --------------------------------------------------
dat_transformed <- transform(dat_raw,
                             transformList(marker_cols, asinh))

## 2.4. Scale from 0 to 1 -------------------------------------------------------

### 2.4.1. Find the minimal and maximal values for each marker
dat_merged <- merge_flowSet(dat_transformed, rename = FALSE)

minimal <- map_dbl(marker_cols, ~ min(dat_merged[, .x]))
names(minimal) <- marker_cols

maximal <- map_dbl(marker_cols, ~ max(dat_merged[, .x]))
names(maximal) <- marker_cols

### 2.5.2. Make a list with transformation functions
scaleTrans <- lapply(marker_cols, function(marker){
  funct <- scaleTransform(transformationId = marker,
                          a = minimal[marker],
                          b = maximal[marker])
  return(funct)
})

### 2.4.2. Scale
dat_scaled <- transform(dat_transformed,
                        transformList(marker_cols, scaleTrans))

## 2.5. Positivity thresholds for markers --------------------------------------

##### The positivity threshold for each marker should be found manually by comparing the samples with unstained sample, FMOs, and/or isotype controls.
##### Here, we only calculate where those previously determined thresholds will be in the transformed and scaled data.

### 2.5.1. Find thresholds
thresholds <- metadata[['panel']]$threshold
names(thresholds) <- marker_cols

thresholds_transformed <- asinh(thresholds)

thresholds_scaled <- map_dbl(marker_cols, function(marker){
  threshold_value <- (thresholds_transformed[marker] - minimal[marker])/(maximal[marker] - minimal[marker])
  return(threshold_value)
})
names(thresholds_scaled) <- marker_cols

### 2.5.2. Plot the positivity thresholds on the histograms (not included in the paper)
pdf(paste0(wd, '010_hist_scaled_gated.pdf'), width = 5, height = 2.5)
map(marker_cols, function(marker){
  p <- ggcyto(dat_scaled[1], aes(x = !!sym(marker))) +
    geom_density() +
    geom_vline(xintercept = thresholds_scaled[marker],
               color = 'blue') +
    ggtitle('After arcsinh transformation and scaling. Vertical line = pos/neg gate') +
    theme_bw()
  print(p)
})
dev.off()

## 2.7. Merge cells from each patient into one object --------------------------
##### It can be done without batch correction since the samples were stained and recorded next to each other
dat_scaled_merged <- merge_flowSet(dat_scaled, rename = FALSE)

## 2.8. Downsample -------------------------------------------------------------
##### We have to do this because hundreds of thousands cells will be too computationally heavy for dimensionality reduction.
##### Take equal amount of cells from each condition
dat_downsampled <- dat_scaled_merged %>% 
  group_by(exp) %>%
  sample_n(size = 5000, replace = FALSE) %>%
  ungroup()

## 2.9. Create a dataframe with only data for PCA/tSNE -------------------------
dat_for_dimred <- dplyr::select(dat_downsampled, -exp)

saveRDS(dat_for_dimred, paste0(wd, '020_dat_for_dimred.rds'))
##### dat_for_dimred <- readRDS(paste0(wd, '020_dat_for_dimred.rds'))

## 2.10. For each marker, make a column that specifies whether a cell is above or below threshold for this marker
for (marker in marker_cols){
  dat_downsampled[, paste0('threshold_', marker)] <- ifelse(
    dat_downsampled[, marker] <= thresholds_scaled[marker],
    'below', 'above'
  )
}

## 2.10. Add columns that specify the patient and LNP condition separately -----
metadata_to_add <- metadata[['samples']] %>%
  dplyr::select(sample_id, condition, patient_id)

dat_downsampled <- dplyr::rename(dat_downsampled, sample_id = exp) %>%
  left_join(metadata_to_add, by = 'sample_id')

# 3. PCA =======================================================================

## 3.1. Perform PCA ------------------------------------------------------------
pca_out <- PCA(
  dat_for_dimred,
  scale.unit = TRUE,
  graph = FALSE,
  ncp = 2
)

saveRDS(pca_out, paste0(wd, '030_pca_res.rds'))
##### pca_out <- readRDS(paste0(wd, '030_pca_res.rds'))

## 3.2. Add PCA coordinates to the data ----------------------------------------
dat_downsampled <- cbind(dat_downsampled, as.data.frame(pca_out$ind$coord)) %>%
  dplyr::rename(PC1 = Dim.1,
                PC2 = Dim.2)

## 4.3. Plot PCA (not included in the paper) -----------------------------------

### 4.3.1. Colored by condition (knockdown or mismatch control)
pdf(paste0(wd, '040_PCA_scaledToVariance_byCond.pdf'), height = 4, width = 6)
ggplot(dat_downsampled %>% sample_frac(),
              aes(x = PC1, y = PC2, color = condition)) +
       geom_point(size = 0.1, alpha = 0.6) +
       theme_bw() +
       scale_color_manual(values = kdmm_palette) +
       labs(x = paste0('PC1, ', as.character(signif(pca_out$eig[1,'percentage of variance'], digits = 3)), '% variance'),
            y = paste0('PC2, ', as.character(signif(pca_out$eig[2, 'percentage of variance'], digits = 3)), '% variance')
       ) +
       theme(axis.text.x = element_blank(),  
             axis.text.y = element_blank(),  
             axis.ticks = element_blank())
pca_list
dev.off()

### 4.3.2. Colored by patient
pdf(paste0(wd, '041_PCA_scaledToVariance_ptsABD_bySample.pdf'), height = 4, width = 4.75)
ggplot(dat_downsampled %>% sample_frac(),
       aes(x = PC1, y = PC2, color = patient_id)) +
  geom_point(size = 0.1, alpha = 0.6) +
  theme_bw() +
  scale_color_manual(values = c(cbPalette[1], cbPalette[3], cbPalette[7])) +
  labs(x = paste0('PC1, ', as.character(signif(pca_out$eig[1,'percentage of variance'], digits = 3)), '% variance'),
       y = paste0('PC2, ', as.character(signif(pca_out$eig[2, 'percentage of variance'], digits = 3)), '% variance')
  ) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())
dev.off()

### 4.3.3. Colored by patient, split by condition
pdf(paste0(wd, '042_PCA_scaledToVariance_ptsABD_colByPt_splitByCond.pdf'), height = 7, width = 4.75)
ggplot(dat_downsampled %>% sample_frac(),
       aes(x = PC1, y = PC2, color = patient_id)) +
  geom_point(size = 0.1, alpha = 0.6) +
  facet_grid(condition ~ .,
             labeller = as_labeller(c(
               'RUNX1::RUNX1T1 knockdown' = 'RUNX1::RUNX1T1 KD',
               'mismatch control' = 'mismatch control'
             ))) +
  theme_bw() +
  scale_color_manual(values = c(cbPalette[1], cbPalette[3], cbPalette[7])) +
  labs(x = paste0('PC1, ', as.character(signif(pca_out$eig[1,'percentage of variance'], digits = 3)), '% variance'),
       y = paste0('PC2, ', as.character(signif(pca_out$eig[2, 'percentage of variance'], digits = 3)), '% variance')
  ) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())
dev.off()

# 4. tSNE ======================================================================

## 4.1. Perform tSNE -----------------------------------------------------------
out_tsne <- Rtsne(as.matrix(dat_for_dimred), 
                  pca = FALSE, 
                  verbose = TRUE)

## 4.2. Add tSNE coordinates to the data ---------------------------------------
dat_downsampled <- cbind(dat_downsampled, as.data.frame(out_tsne$Y)) %>%
  dplyr::rename(tSNE1 = V1,
                tSNE2 = V2)

saveRDS(dat_downsampled, paste0(wd, '990_data_forPlotting.rds'))
##### dat_downsampled <- readRDS(paste0(wd, '990_data_forPlotting.rds'))

## 4.3. Plot tSNE -------------------------------------------------------------- 

### 4.3.1. Colored by patient, split by condition (Suppl. Figure 6C)
pdf(paste0(wd, '047_tSNE_ptABD_colByPt_splitByCond.pdf'), height = 5.5, width = 4)
ggplot(dat_downsampled %>% sample_frac(),
       aes(x = tSNE1, y = tSNE2, color = patient_id)) +
  geom_point(size = 0.1) +
  facet_grid(condition ~ .,
             labeller = as_labeller(c(
               'RUNX1::RUNX1T1 knockdown' = 'RUNX1::RUNX1T1 KD',
               'mismatch control' = 'mismatch control'
             ))) +
  theme_bw() +
  scale_color_manual(values = c(cbPalette[1], cbPalette[3], cbPalette[7])) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())
dev.off()

### 4.3.4. Colored by expression of markers, split by condition (Suppl. Figure 6B)
pdf(paste0(wd, '051_tSNE_byMarkers.pdf'),
    height = 10, width = 10)

plotlist <- map(marker_cols, function(marker){
  p <- ggplot(dat_downsampled,
              aes(x = tSNE1, y = tSNE2)) +
    geom_point(size = 0.1,
               aes(colour = !!sym(marker), group = paste0('threshold_', marker)),
               data = subset(dat_downsampled, get(marker) > thresholds_scaled[marker])) +
    geom_point(size = 0.1,
               color = 'grey',
               data = subset(dat_downsampled, get(marker) <= thresholds_scaled[marker])) +
    scale_color_viridis_c(begin = 0,
                          end = 1,
                          direction = 1,
                          na.value = 'grey',
                          name = paste0(marker, '\nscaled expression')) +
    facet_grid(condition ~ .,
               labeller = as_labeller(c(
                 'RUNX1::RUNX1T1 knockdown' = 'RUNX1::RUNX1T1 KD',
                 'mismatch control' = 'mismatch control'
               ))) +
    theme_bw() +
    theme(legend.title = element_text(size = 8),
          legend.position = 'right',
          axis.text.x = element_blank(),  
          axis.text.y = element_blank(),  
          axis.ticks = element_blank())
  return(p)
})
names(plotlist) <- marker_cols
myplot <- ggarrange(plotlist = plotlist, ncol = 3, nrow = 3)
print(myplot)
dev.off()

# 6. Plot the 'double exposure' tSNEs for CD34 and CD38 (Suppl. Figure 6D) =====
dat_downsampled$CD34_adj <- dat_downsampled[, 'CD34-APC'] - thresholds_scaled['CD34-APC']
dat_downsampled$CD34_adj[dat_downsampled$CD34_adj < 0] <- 0
max_cd34 <- max(dat_downsampled$CD34_adj)
dat_downsampled$CD34_adj <- dat_downsampled$CD34_adj/max_cd34

dat_downsampled$CD38_adj <- dat_downsampled[, 'CD38-PE_Dazzle594'] - thresholds_scaled['CD38-PE_Dazzle594']
dat_downsampled$CD38_adj[dat_downsampled$CD38_adj < 0] <- 0
max_cd38 <- max(dat_downsampled$CD38_adj)
dat_downsampled$CD38_adj <- dat_downsampled$CD38_adj/max_cd38

dat_downsampled$cd34_38_color <- mapply(function(x,y) rgb(0, x, y), dat_downsampled$CD34_adj, dat_downsampled$CD38_adj)

saveRDS(dat_downsampled, paste0(wd, '990_data_forPlotting.rds'))
##### dat_downsampled <- readRDS(paste0(wd, '990_data_forPlotting.rds'))

pdf(paste0(wd, '070_tSNE_CD34-38_doubleExposure.pdf'),
    height = 5.5, width = 6)

##### manual legend
legend_space <- expand.grid(CD38 = seq(0, 1, length.out = 10),
                           CD34 = seq(0, 1, length.out = 10))
legend_space$color <- apply(legend_space, 1, function(x){
  rgb(0, x[2], x[1], maxColorValue = 1)
})

legend_plot <- ggplot(legend_space,
                      aes(x = CD38, y = CD34, fill = I(color))) +
  geom_tile() +
  coord_fixed() +
  theme_bw() +
  labs(x = 'CD38 scaled expression',
       y = 'CD34 scaled expression')

tsne <- ggplot(dat_downsampled,
               aes(x = tSNE1, y = tSNE2, color = I(cd34_38_color))) +
  geom_point(size = 0.1) +
  facet_grid(condition ~ .,
             labeller = as_labeller(c(
               'RUNX1::RUNX1T1 knockdown' = 'RUNX1::RUNX1T1 KD',
               'mismatch control' = 'mismatch control'
             ))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())

ggarrange(tsne, legend_plot, nrow = 1, ncol = 2)
  
dev.off()

# 99. Session info =============================================================
sink(paste0(wd, '999_sessionInfo.txt'))
sessionInfo()
sink()