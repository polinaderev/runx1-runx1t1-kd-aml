library(statmod)
library(tidyverse)

palette <- c('#222222', '#2538a5', '#f41626')

set.seed(42)

##### directory for output files, replace with your path
wd <- 'C:/polina/analysis/pderevianko/runx1-runx1t1-kd-aml/ELDA/out/'

# 1. Submit the data ===========================================================

##### number of positive cases out of tested trials
response <- c(6,6,6,6,6,5,6,6,3,5,5,3,5,4,1)

##### expected number of cells in assay
dose <- c(4000,4000,4000,2000,2000,2000,1000,1000,1000,500,500,500,250,250,250)

##### number of trials at each dose
tested <- rep(6, length(response))

##### to which experimental condition the responses belong
group <- rep(c('Mock', 'siMM', 'siRE'), length(unique(dose)))

# 2. Perform the ELDA ==========================================================
res <- elda(
  response = response,
  dose = dose,
  tested = tested,
  group = group,
  observed = TRUE,
  confidence = 0.95,
  test.unit.slope = TRUE
)

plot(res)

# 3. Plot ======================================================================

## 3.1. Make the necessary dataframes for ggplot2 ------------------------------
df <- data.frame(
  dose = res$dose,
  response = res$response,
  tested = res$tested,
  group = factor(res$group)
) %>%
  mutate(
    full_response = response == tested,
    adj_response = ifelse(full_response, response-0.5, response),
    log_nonres = log(1 - adj_response/tested),  
    pch_type = ifelse(full_response, "100%", "<100%")
  )

ci_df <- data.frame(
  group = levels(df$group),
  CI_lower = res$CI[, 1],
  CI_median = res$CI[, 2],
  CI_upper = res$CI[, 3]
)

## 3.2. Plot -------------------------------------------------------------------
p <- ggplot(df, 
            aes(x = dose, y = log_nonres, color = group)) +
  geom_point(aes(shape = pch_type), size = 3) +
  scale_shape_manual(values = c("100%" = 1, "<100%" = 6)) +
  labs(
    x = "Number of cells seeded",
    y = "ln(fraction of wells without cells)",
    shape = "Percentage of wells\nwith cell growth",
    color = "Group"
  ) +
  theme_minimal() +
  scale_color_manual(values = palette)

for (i in 1:nrow(ci_df)) {
  p <- p +
    geom_abline(slope = -1/ci_df$CI_median[i], intercept = 0,
                color = palette[i], size = 0.5) +
    geom_abline(slope = -1/ci_df$CI_lower[i], intercept = 0,
                linetype = "dashed", color = palette[i], linewidth = 0.4) +
    geom_abline(slope = -1/ci_df$CI_upper[i], intercept = 0,
                linetype = "dashed", color = palette[i], linewidth = 0.4)
}

pdf(paste0(wd, '010_ELDAplot.pdf'), width = 6, height = 4)
p
dev.off()
