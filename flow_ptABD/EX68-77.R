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

in_path_68 <- 'flow/in/01_flowjo_EX68/' ##### Experiment with the 3 patient samples
in_path_77 <- 'flow/in/01_flowjo_EX77/' ##### Experiment with AML PDX
in_path <- c(in_path_68, in_path_77)
names(in_path) <- c('EX68', 'EX77')
wd <- 'flow/out/' ##### output directory

##### reproducibility
set.seed(42)

# 1. Load in the data ==========================================================

## 1.1. Metadata Excel file that contains info about samples and about the panel -----
excel_path <- map_chr(in_path, ~ paste0(.x, 'metadata.xlsx'))
metadata_sheets <- map(excel_path, excel_sheets)

metadata <- lapply(names(metadata_sheets), function(experiment_name){
  list_new <- map(metadata_sheets[[experiment_name]], 
                  ~ read_excel(excel_path[[experiment_name]], 
                               sheet = .x))
  names(list_new) <- metadata_sheets[[experiment_name]]
  return(list_new)
})
names(metadata) <- names(in_path)

## 1.2. FCS files --------------------------------------------------------------

### 1.2.1. Load 
dat_raw <- lapply(names(in_path), function(experiment_name){
  flowset <- read.flowSet(files = metadata[[experiment_name]][['samples']]$file_name,
                          path = in_path[experiment_name],
                          descriptions = metadata[[experiment_name]][['samples']]$sample_id,
                          transformation = FALSE,
                          truncate_max_range = FALSE)
  sampleNames(flowset) <- metadata[[experiment_name]][['samples']]$sample_id
  pData(flowset)$name <- metadata[[experiment_name]][['samples']]$sample_id
  return(flowset)
})
names(dat_raw) <- names(in_path)

### 1.2.2. Replace the column names (channels) with marker names
rename_vector <- lapply(names(in_path), function(experiment_name){
  vec <- setNames(metadata[[experiment_name]][['panel']]$marker, 
                  metadata[[experiment_name]][['panel']]$channel)
  return(vec)
})
names(rename_vector) <- names(in_path)

columns_to_rename <- lapply(names(in_path), function(experiment_name){
  vec <- colnames(dat_raw[[experiment_name]]) %in% names(rename_vector[[experiment_name]])
  return(vec)
})
names(columns_to_rename) <- names(in_path)

dat_raw <- lapply(names(in_path), function(experiment_name){
  colnames(dat_raw[[experiment_name]])[columns_to_rename[[experiment_name]]] <- rename_vector[[experiment_name]][colnames(dat_raw[[experiment_name]])[columns_to_rename[[experiment_name]]]]
  return(dat_raw[[experiment_name]])
})
names(dat_raw) <- names(in_path)

# 2. Prepare the data ==========================================================

## 2.1. Only protein marker columns will be used for the analysis --------------
dat_raw <- lapply(names(in_path), function(experiment_name){
  marker_cols <- metadata[[experiment_name]][['panel']]$marker
  dat_raw[[experiment_name]] <- dat_raw[[experiment_name]][, marker_cols]
  return(dat_raw[[experiment_name]])
})
names(dat_raw) <- names(in_path)

marker_cols <- lapply(metadata, function(sublist){
  vec <- sublist[['panel']]$marker
  return(vec)
})

## 2.2. Remove cells with negative values --------------------------------------

### 2.2.1. Cells per sample before filtering
fsApply(dat_raw[[1]], nrow)
# [,1]
# A_GMCSF_RE 49974
# A_GMCSF_MM 41623
# B_GMCSF_RE 32246
# B_GMCSF_MM 47754
# D_GMCSF_RE 53371
# D_GMCSF_MM 42701

fsApply(dat_raw[[2]], nrow)
# [,1]
# C_GMCSF_RE 29437
# C_GMCSF_MM 50787

### 2.2.2. Prepare gates for each marker of non-negative values
gates <- lapply(marker_cols, function(vec){
  vec_new <- map(vec,
                     function(marker){
                       mat <- matrix(c(0,Inf), ncol = 1, dimnames = list(c('min', 'max'), marker))
                       gate <- rectangleGate(filterId = paste0(sub('-.*', '', marker), '_clean'), .gate = mat)
                       return(gate)
                       })
  return(vec_new)
})

### 2.2.3. Gate and subset 
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

### 2.2.4. Cells per sample after filtering
fsApply(dat_filt[[1]], nrow)
# [,1]
# A_GMCSF_RE 43055
# A_GMCSF_MM 36663
# B_GMCSF_RE 27294
# B_GMCSF_MM 40580
# D_GMCSF_RE 45190
# D_GMCSF_MM 36145

fsApply(dat_filt[[2]], nrow)
# [,1]
# C_GMCSF_RE 27982
# C_GMCSF_MM 48776

## 2.3. Arcsinh transformation --------------------------------------------------

### 2.3.1. Apply transformation
dat_transformed <- lapply(names(in_path), function(experiment_name){
  dat_transformed <- transform(dat_filt[[experiment_name]],
                               transformList(marker_cols[[experiment_name]], asinh))
  return(dat_transformed)
})
names(dat_transformed) <- names(in_path)

## 2.4. Subset to samples that will be taken together for plotting -------------

##### In this particular analysis, I'm only interested in samples with GMCSF.
##### By patient and all primary samples together.
ptA <- sampleNames(dat_transformed[['EX68']])[grepl('^A_GMCSF', sampleNames(dat_transformed[['EX68']]))]
ptB <- sampleNames(dat_transformed[['EX68']])[grepl('^B_GMCSF', sampleNames(dat_transformed[['EX68']]))]
ptC <- sampleNames(dat_transformed[['EX77']])[grepl('^C_GMCSF', sampleNames(dat_transformed[['EX77']]))]
ptD <- sampleNames(dat_transformed[['EX68']])[grepl('^D_GMCSF', sampleNames(dat_transformed[['EX68']]))]
ptABD <- sampleNames(dat_transformed[['EX68']])

samples_for_analysis <- list('patientA' = ptA, 
                             'patientB' = ptB, 
                             'patientC_PDX' = ptC,
                             'patientD' = ptD,
                             'patientsABD' = ptABD)

sample_names <- c('patientA', 'patientB', 'patientC_PDX', 'patientD', 'patientsABD')

dat_subset <- lapply(samples_for_analysis, function(sampleIDs){
  if (all(sampleIDs %in% sampleNames(dat_transformed[['EX68']]))){
    dat_new <- dat_transformed[['EX68']][sampleIDs]
  } else {
    dat_new <- dat_transformed[['EX77']][sampleIDs]
  }
  return(dat_new)
})

## 2.5. Scale from 0 to 1 -------------------------------------------------------

##### Since the flow cytometry for patient C was recorded separately on a different day than the other 3 patients, I cannot scale them to the same maximal and minimal value.

### 2.5.1. Find the minimal and maximal values for each marker
dat_merged <- lapply(dat_subset, function(flowset){
  merged <- merge_flowSet(flowset, rename = FALSE)
  return(merged)
})

minimal <- lapply(dat_merged, function(df){
  minimal <- map_dbl(marker_cols[[1]], ~ min(df[, .x]))
  names(minimal) <- marker_cols[[1]]
  return(minimal)
})

maximal <- lapply(dat_merged, function(df){
  maximal <- map_dbl(marker_cols[[1]], ~ max(df[, .x]))
  names(maximal) <- marker_cols[[1]]
  return(maximal)
})

### 2.5.2. Make a list with transformation functions
scaleTrans <- lapply(names(samples_for_analysis), function(patient){
  scaleTrans <- lapply(marker_cols[[1]], function(marker){
    funct <- scaleTransform(transformationId = marker,
                            a = minimal[[patient]][marker],
                            b = maximal[[patient]][marker])
    return(funct)
  })
  names(scaleTrans) <- marker_cols[[1]]
  return(scaleTrans)
})
names(scaleTrans) <- names(samples_for_analysis)

### 2.5.3. Scale
dat_scaled <- lapply(names(samples_for_analysis), function(patient){
  scaled <- transform(dat_subset[[patient]],
                      transformList(marker_cols[[1]], scaleTrans[[patient]]))
  return(scaled)
})
names(dat_scaled) <- names(samples_for_analysis)

## 2.6. Positivity thresholds for markers --------------------------------------

##### The positivity threshold for each marker should be found manually by comparing the samples with unstained sample, FMOs, and/or isotype controls.
##### Here, we only calculate where those previously determined thresholds will be in the transformed and scaled data.

### 2.5.1. Find thresholds
thresholds_68 <- metadata[['EX68']][['panel']]$threshold
names(thresholds_68) <- metadata[['EX68']][['panel']]$marker

thresholds_77 <- metadata[['EX77']][['panel']]$threshold
names(thresholds_77) <- metadata[['EX77']][['panel']]$marker

thresholds <- list(thresholds_68, thresholds_68, thresholds_77, thresholds_68, thresholds_68)
names(thresholds) <- names(samples_for_analysis)

thresholds_transformed <- map(thresholds, asinh)

thresholds_scaled <- lapply(names(samples_for_analysis), function(patient){
  scaled <- map_dbl(marker_cols[[1]], function(marker){
    threshold_value <- (thresholds_transformed[[patient]][marker] - minimal[[patient]][marker])/(maximal[[patient]][marker] - minimal[[patient]][marker])
    return(threshold_value)
  })
  names(scaled) <- marker_cols[[1]]
  return(scaled)
})
names(thresholds_scaled) <- names(samples_for_analysis)

### 2.5.2. Plot the positivity thresholds on the histograms
for (patient in names(dat_scaled)){
  pdf(paste0(wd, '010_hist_scaled_gated_', patient, '.pdf'), width = 5, height = 2.5)
  map(marker_cols[[1]], function(marker){
    p <- ggcyto(dat_scaled[[patient]], aes(x = !!sym(marker))) +
      geom_density() +
      geom_vline(xintercept = thresholds_scaled[[patient]][marker],
                 color = 'blue') +
      ggtitle('After filtering of non-negative values, arcsinh transformation and scaling. Vertical line = pos/neg gate') +
      theme_bw()
    print(p)
  })
  dev.off()
}

## 2.7. Merge cells from each patient into one object --------------------------
##### It can be done without batch correction since the samples were stained and recorded next to each other
dat_scaled_merged <- map(dat_scaled, ~ merge_flowSet(.x, rename = FALSE))

## 2.8. Downsample -------------------------------------------------------------
##### We have to do this because hundreds of thousands cells will be too computationally heavy for dimensionality reduction.
##### Take equal amount of cells from each condition
dat_downsampled <- map(dat_scaled_merged, 
                         ~ group_by(.x, exp) %>%
                         sample_n(size = 5000, replace = FALSE) %>%
                         ungroup())

## 2.9. Create a dataframe with only data for PCA/tSNE -------------------------
dat_for_dimred <- map(dat_downsampled, 
                    ~ dplyr::select(.x, -exp))

saveRDS(dat_for_dimred, paste0(wd, '020_dat_for_dimred.rds'))
##### dat_for_dimred <- readRDS(paste0(wd, '020_dat_for_dimred.rds'))

## 2.10. For each marker, make a column that specifies whether a cell is above or below threshold for this marker
dat_downsampled <- lapply(names(dat_downsampled), function(patient){
  for (marker in marker_cols[[1]]){
    dat_downsampled[[patient]][, paste0('threshold_', marker)] <- ifelse(
      dat_downsampled[[patient]][, marker] <= thresholds_scaled[[patient]][marker],
      'below', 'above'
    )
  }
  return(dat_downsampled[[patient]])
})
names(dat_downsampled) <- names(samples_for_analysis)
# 
## 2.11. Add columns that specify the patient and LNP condition separately -----
metadata_to_add_68 <- metadata[['EX68']][['samples']] %>%
  dplyr::select(sample_id, condition, patient_id)

metadata_to_add_77 <- metadata[['EX77']][['samples']] %>%
  dplyr::select(sample_id, condition, patient_id)

metadata_to_add <- list(metadata_to_add_68, metadata_to_add_68, metadata_to_add_77, metadata_to_add_68, metadata_to_add_68)
names(metadata_to_add) <- names(samples_for_analysis)

dat_downsampled <- lapply(names(samples_for_analysis), function(patient){
  df <- dplyr::rename(dat_downsampled[[patient]], sample_id = exp) %>%
    left_join(metadata_to_add[[patient]], by = 'sample_id')
    return(df)
})
names(dat_downsampled) <- names(samples_for_analysis)

# 3. PCA =======================================================================

## 3.1. Perform PCA ------------------------------------------------------------
pca_out <- map(dat_for_dimred,
               ~ PCA(
                 .x,
                 scale.unit = TRUE,
                 graph = FALSE,
                 ncp = 2
               ))

saveRDS(pca_out, paste0(wd, '030_pca_res.rds'))
##### pca_out <- readRDS(paste0(wd, '030_pca_res.rds'))

## 3.2. Add PCA coordinates to the data ----------------------------------------
dat_downsampled <- map2(dat_downsampled, pca_out,
                        ~ cbind(.x, as.data.frame(.y$ind$coord)) %>%
                          dplyr::rename(PC1 = Dim.1,
                                        PC2 = Dim.2))

## 4.3. Plot PCA ---------------------------------------------------------------

### 4.3.1. Colored by condition (knockdown or mismatch control)
pdf(paste0(wd, '040_PCA_scaledToVariance_byCond.pdf'), height = 4, width = 6)
pca_list <- pmap(list(dat_downsampled, names(dat_downsampled), pca_out),
     ~ ggplot(..1 %>% sample_frac(),
              aes(x = PC1, y = PC2, color = condition)) +
       geom_point(size = 0.1, alpha = 0.6) +
       theme_bw() +
       scale_color_manual(values = kdmm_palette) +
       ggtitle(paste0(..2)) +
       labs(x = paste0('PC1, ', as.character(signif(..3$eig[1,'percentage of variance'], digits = 3)), '% variance'),
            y = paste0('PC2, ', as.character(signif(..3$eig[2, 'percentage of variance'], digits = 3)), '% variance')
       ) +
       theme(axis.text.x = element_blank(),  
             axis.text.y = element_blank(),  
             axis.ticks = element_blank())
)
pca_list
dev.off()

### 4.3.2. Colored by patient
pdf(paste0(wd, '041_PCA_scaledToVariance_ptsABD_bySample.pdf'), height = 4, width = 4.75)
ggplot(dat_downsampled[['patientsABD']] %>% sample_frac(),
       aes(x = PC1, y = PC2, color = patient_id)) +
  geom_point(size = 0.1, alpha = 0.6) +
  theme_bw() +
  scale_color_manual(values = c(cbPalette[1], cbPalette[3], cbPalette[7])) +
  ggtitle(names(dat_downsampled)[5]) +
  labs(x = paste0('PC1, ', as.character(signif(pca_out[['patientsABD']]$eig[1,'percentage of variance'], digits = 3)), '% variance'),
       y = paste0('PC2, ', as.character(signif(pca_out[['patientsABD']]$eig[2, 'percentage of variance'], digits = 3)), '% variance')
  ) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())
dev.off()

### 4.3.3. Colored by patient, split by condition
pdf(paste0(wd, '042_PCA_scaledToVariance_ptsABD_colByPt_splitByCond.pdf'), height = 7, width = 4.75)
ggplot(dat_downsampled[['patientsABD']] %>% sample_frac(),
       aes(x = PC1, y = PC2, color = patient_id)) +
  geom_point(size = 0.1, alpha = 0.6) +
  facet_grid(condition ~ .,
             labeller = as_labeller(c(
               'RUNX1::RUNX1T1 knockdown' = 'RUNX1::RUNX1T1 KD',
               'mismatch control' = 'mismatch control'
             ))) +
  theme_bw() +
  scale_color_manual(values = c(cbPalette[1], cbPalette[3], cbPalette[7])) +
  ggtitle(names(dat_downsampled)[5]) +
  labs(x = paste0('PC1, ', as.character(signif(pca_out[['patientsABD']]$eig[1,'percentage of variance'], digits = 3)), '% variance'),
       y = paste0('PC2, ', as.character(signif(pca_out[['patientsABD']]$eig[2, 'percentage of variance'], digits = 3)), '% variance')
  ) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())
dev.off()

# 4. tSNE ======================================================================

## 4.1. Perform tSNE -----------------------------------------------------------
out_tsne <- map(dat_for_dimred,
                ~ Rtsne(as.matrix(.x), 
                  pca = FALSE, 
                  verbose = TRUE))

## 4.2. Add tSNE coordinates to the data ---------------------------------------
dat_downsampled <- map2(dat_downsampled, out_tsne,
                        ~ cbind(.x, as.data.frame(.y$Y)) %>%
                          dplyr::rename(tSNE1 = V1,
                                        tSNE2 = V2))

saveRDS(dat_downsampled, paste0(wd, '990_data_forPlotting.rds'))
##### dat_downsampled <- readRDS(paste0(wd, '990_data_forPlotting.rds'))

## 4.3. Plot tSNE -------------------------------------------------------------- 

### 4.3.1. Colored by condition
pdf(paste0(wd, '045_tSNE_byCond.pdf'), height = 4, width = 6)
tsne_list <- map2(dat_downsampled, names(dat_downsampled),
                 ~ ggplot(.x %>% sample_frac(),
                          aes(x = tSNE1, y = tSNE2, color = condition)) +
                   geom_point(size = 0.15) +
                   theme_bw() +
                   scale_color_manual(values = kdmm_palette) +
                   ggtitle(paste0(.y)) +
                   theme(axis.text.x = element_blank(),  
                         axis.text.y = element_blank(),  
                         axis.ticks = element_blank())
)
tsne_list
dev.off()

### 4.3.2. Colored by patient
pdf(paste0(wd, '046_tSNE_ptABD_byPt.pdf'), height = 4, width = 4.5)
ggplot(dat_downsampled[['patientsABD']] %>% sample_frac(),
       aes(x = tSNE1, y = tSNE2, color = patient_id)) +
  geom_point(size = 0.15) +
  theme_bw() +
  scale_color_manual(values = c(cbPalette[1], cbPalette[3], cbPalette[7])) +
  ggtitle(names(dat_downsampled)[4]) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())
dev.off()

### 4.3.3. Colored by patient, split by condition
pdf(paste0(wd, '047_tSNE_ptABD_colByPt_splitByCond.pdf'), height = 5.5, width = 4)
ggplot(dat_downsampled[['patientsABD']] %>% sample_frac(),
       aes(x = tSNE1, y = tSNE2, color = patient_id)) +
  geom_point(size = 0.1) +
  facet_grid(condition ~ .,
             labeller = as_labeller(c(
               'RUNX1::RUNX1T1 knockdown' = 'RUNX1::RUNX1T1 KD',
               'mismatch control' = 'mismatch control'
             ))) +
  theme_bw() +
  scale_color_manual(values = c(cbPalette[1], cbPalette[3], cbPalette[7])) +
  ggtitle(names(dat_downsampled)[5]) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(),  
        axis.ticks = element_blank())
dev.off()

### 4.3.4. Colored by expression of markers, split by condition
tsne_list <- vector('list', length(samples_for_analysis))
names(tsne_list) <- names(samples_for_analysis)

for (patient in names(samples_for_analysis)){
  
  pdf(paste0(wd, '050_tSNE_byMarkers_', patient, '.pdf'),
      height = 10, width = 10)
  
  plotlist <- map(marker_cols[[1]], function(marker){
    p <- ggplot(dat_downsampled[[patient]],
                aes(x = tSNE1, y = tSNE2)) +
      geom_point(size = 0.1,
                 aes(colour = !!sym(marker), group = paste0('threshold_', marker)),
                 data = subset(dat_downsampled[[patient]], get(marker) > thresholds_scaled[[patient]][marker])) +
      geom_point(size = 0.1,
                 color = 'grey',
                 data = subset(dat_downsampled[[patient]], get(marker) <= thresholds_scaled[[patient]][marker])) +
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
  names(plotlist) <- marker_cols[[1]]
  tsne_list[[patient]] <- plotlist
  myplot <- ggarrange(plotlist = plotlist, ncol = 3, nrow = 3)
  print(myplot)
  dev.off()
}

# 5. Density plots for markers =================================================
density_list <- vector('list', length(samples_for_analysis))
names(density_list) <- names(samples_for_analysis)
  
for (patient in names(samples_for_analysis)){
  
  pdf(paste0(wd, '060_density_byMarkers_', patient, '.pdf'),
      height = 4, width = 10)
  
  plotlist <- map(marker_cols[[1]], function(marker){
    p <- ggplot(dat_downsampled[[patient]],
                aes(x = !!sym(marker), fill = condition, color = condition)) +
      geom_density(alpha = 0.4, size = 0.5) +
      theme_bw() +
      theme(legend.position = 'none') +
      scale_fill_manual(values = kdmm_palette) +
      scale_color_manual(values = kdmm_palette) +
      xlab(paste0(marker, ' scaled intensity'))
  })
  names(plotlist) <- marker_cols[[1]]
  density_list[[patient]] <- plotlist
  myplot <- ggarrange(plotlist = plotlist, ncol = 3, nrow = 3)
  print(myplot)
  dev.off()
}

# 6. Plot the 'double exposure' tSNEs for CD34 and CD38 ========================
dat_downsampled <- lapply(names(samples_for_analysis), function(patient){
  df <- dat_downsampled[[patient]]
  
  df$CD34_adj <- df[, 'CD34-APC'] - thresholds_scaled[[patient]]['CD34-APC']
  df$CD34_adj[df$CD34_adj < 0] <- 0
  max_cd34 <- max(df$CD34_adj)
  df$CD34_adj <- df$CD34_adj/max_cd34
  
  df$CD38_adj <- df[, 'CD38-PE_Dazzle594'] - thresholds_scaled[[patient]]['CD38-PE_Dazzle594']
  df$CD38_adj[df$CD38_adj < 0] <- 0
  max_cd38 <- max(df$CD38_adj)
  df$CD38_adj <- df$CD38_adj/max_cd38
  
  df$cd34_38_color <- mapply(function(x,y) rgb(0, x, y), df$CD34_adj, df$CD38_adj)
  return(df)
})
names(dat_downsampled) <- names(samples_for_analysis)

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

tsne_de_list <- lapply(sample_names, function(patient){
  tsne <- ggplot(dat_downsampled[[patient]],
              aes(x = tSNE1, y = tSNE2, color = I(cd34_38_color))) +
    geom_point(size = 0.1) +
    facet_grid(condition ~ .,
               labeller = as_labeller(c(
                 'RUNX1::RUNX1T1 knockdown' = 'RUNX1::RUNX1T1 KD',
                 'mismatch control' = 'mismatch control'
               ))) +
    theme_bw() +
    ggtitle(patient) +
    theme(axis.text.x = element_blank(),  
          axis.text.y = element_blank(),  
          axis.ticks = element_blank())
  
  p <- ggarrange(tsne, legend_plot, nrow = 1, ncol = 2)
  return(p)
})

names(tsne_de_list) <- sample_names

print(tsne_de_list)

dev.off()

# 99. Session info =============================================================
sink(paste0(wd, '999_sessionInfo.txt'))
sessionInfo()
sink()