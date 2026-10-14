# Unknown content (as Fig. 5C of the omnideconv benchmark) with the results of the rectangle revision:
# lambrechts split by donor, tumour cells as unknown content, eight omnideconv methods, rectangle and SECRET.
# Based on fig5/fig_5C.R; the computation of the metrics is the same.
#
# Usage: Rscript fig_unknown_content.R [<results_dir>] [<output_dir>]

library(tidyverse)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else '/nfs/data/omnideconv_benchmarking_clean/benchmark_results/rectangle_revision/unknown_content'
output_dir <- if (length(args) >= 2) args[2] else '/nfs/proj/omnideconv_benchmarking/omnideconv/rectangle_revision/figures'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_names <- c(autogenes = 'AutoGeneS', bayesprism = 'BayesPrism', bisque = 'Bisque', cibersortx = 'CIBERSORTx',
                  dwls = 'DWLS', music = 'MuSiC', rectangle = 'Rectangle', scaden = 'Scaden', scdc = 'SCDC', secret = 'SECRET')
vector.frac <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9)

# 1: List result files, methods, replicates and fractions

files <- list.files(results_dir, full.names = FALSE, recursive = TRUE)
files <- files[grepl('replicate', files) & grepl('lambrechts', files) & grepl('\\.rds$', files)]

metadata.table <- tibble(path = files,
                         method_dataset = map_chr(files, function(x) strsplit(x, split = '/')[[1]][1]),
                         replicate = map_chr(files, function(x) strsplit(x, split = '/')[[1]][2]),
                         unknown_fraction = map_chr(files, function(x) strsplit(x, split = '/')[[1]][3])) %>%
  mutate(method_dataset = gsub('_unknown_content_sim', '', method_dataset),
         replicate = as.numeric(gsub('replicate_', '', replicate)),
         unknown_fraction = as.numeric(gsub('deconvolution_Tumor cells_|.rds', '', unknown_fraction))) %>%
  separate(method_dataset, c('method', 'dataset'), sep = '_') %>%
  filter(unknown_fraction %in% vector.frac)

cat('result files per method:\n'); print(table(metadata.table$method))

# 2: Estimated and true fractions of the known cell types in one table. SECRET also estimates an
# 'unknown' proportion; it is kept apart and compared with the true tumour fraction.

read_one <- function(i) {
  res <- readRDS(file.path(results_dir, metadata.table$path[i]))
  info <- metadata.table[i, c('method', 'dataset', 'replicate', 'unknown_fraction')]

  estimated <- as.data.frame(res$deconvolution) %>%
    tibble::rownames_to_column(var = 'sample') %>%
    gather(key = 'celltype', value = 'predicted_value', -sample) %>%
    # Scaden returns the cell-type names with dots instead of spaces
    mutate(celltype = gsub('\\.', ' ', gsub('xxxx', ' ', celltype)))

  truth <- as.data.frame(res$true_cell_fractions)
  tumor <- tibble(sample = colnames(truth), true_tumor = as.numeric(truth['Tumor cells', ]))
  truth <- truth[row.names(truth) != 'Tumor cells', ] %>%
    tibble::rownames_to_column(var = 'celltype') %>%
    gather(key = 'sample', value = 'true_value', -celltype) %>%
    mutate(celltype = gsub('xxxx', ' ', celltype))

  known <- inner_join(estimated, truth, by = c('sample', 'celltype')) %>% bind_cols(info)
  unknown <- estimated %>% filter(celltype == 'unknown') %>% inner_join(tumor, by = 'sample') %>% bind_cols(info)
  list(known = known, unknown = unknown)
}
all_results <- lapply(seq_len(nrow(metadata.table)), read_one)
data <- bind_rows(lapply(all_results, `[[`, 'known'))
secret_unknown <- bind_rows(lapply(all_results, `[[`, 'unknown'))

# 3: Correlation and RMSE per cell type and over all cell types, for every replicate

corelation.results <- data %>%
  group_by(celltype, method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(corr = suppressWarnings(cor(true_value, predicted_value, method = 'pearson')), .groups = 'drop')
corelation.results.all <- data %>%
  group_by(method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(corr = suppressWarnings(cor(true_value, predicted_value, method = 'pearson')), .groups = 'drop') %>%
  mutate(celltype = 'all-cells')
rmse.results <- data %>%
  group_by(celltype, method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(rmse = sqrt(mean((true_value - predicted_value)^2)), .groups = 'drop')
rmse.results.all <- data %>%
  group_by(method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(rmse = sqrt(mean((true_value - predicted_value)^2)), .groups = 'drop') %>%
  mutate(celltype = 'all-cells')

corelation.results <- bind_rows(corelation.results, corelation.results.all)
rmse.results <- bind_rows(rmse.results, rmse.results.all)

# median and interquartile range over the replicates
summarise_replicates <- function(df, value) {
  df %>%
    group_by(celltype, method, unknown_fraction, dataset) %>%
    dplyr::summarize(mean = mean(.data[[value]], na.rm = TRUE),
                     median = median(.data[[value]], na.rm = TRUE),
                     q25 = quantile(.data[[value]], 0.25, na.rm = TRUE),
                     q75 = quantile(.data[[value]], 0.75, na.rm = TRUE),
                     sd = sd(.data[[value]], na.rm = TRUE), .groups = 'drop') %>%
    mutate(method_name = factor(method_names[method], levels = method_names))
}
correlation.results.summary <- summarise_replicates(corelation.results, 'corr')
rmse.results.summary <- summarise_replicates(rmse.results, 'rmse')

write_csv(correlation.results.summary, file.path(output_dir, 'unknown_content_correlation.csv'))
write_csv(rmse.results.summary, file.path(output_dir, 'unknown_content_rmse.csv'))

# 4: Plots

cell_palette <- c('B cells' = '#999933',
                  'Macrophages' = '#CC6677',
                  'T cells CD4 conv' = '#117733',
                  'Stromal cells' = '#288BA8',
                  'all-cells' = '#B3B3B3')

trend_plot <- function(summary, ylab, title, reference_lines) {
  ggplot(summary, aes(x = as.character(unknown_fraction), y = median, group = celltype)) +
    geom_hline(yintercept = reference_lines, linetype = 'dashed', colour = 'grey60') +
    geom_vline(xintercept = '0', linetype = 'dashed', colour = 'grey60') +
    geom_ribbon(aes(fill = celltype, ymin = q25, ymax = q75), alpha = .2) +
    geom_errorbar(aes(ymin = q25, ymax = q75), width = .1) +
    geom_line(aes(color = celltype)) +
    geom_point(aes(color = celltype)) +
    facet_wrap(~method_name, ncol = 5) +
    theme_minimal() +
    theme(legend.position = 'bottom', legend.title = element_blank()) +
    scale_color_manual(values = cell_palette) +
    scale_fill_manual(values = cell_palette) +
    xlab('Unknown cell fraction') +
    ylab(ylab) +
    labs(title = title)
}

corr.ribbon <- trend_plot(correlation.results.summary, 'Pearson correlation',
                          'Correlation trend with increasing tumor fraction', c(0, 0.5))
rmse.ribbon <- trend_plot(rmse.results.summary, 'RMSE',
                          'RMSE trend with increasing tumor fraction', 0)
ggsave(file.path(output_dir, 'unknown_content_correlation.pdf'), corr.ribbon, width = 16, height = 7.5)
ggsave(file.path(output_dir, 'unknown_content_rmse.pdf'), rmse.ribbon, width = 16, height = 7.5)

# SECRET: estimated unknown proportion against the true tumour fraction
if (nrow(secret_unknown) > 0) {
  secret.plot <- ggplot(secret_unknown, aes(x = true_tumor, y = predicted_value)) +
    geom_abline(linetype = 'dashed', colour = 'grey60') +
    geom_jitter(width = 0.008, height = 0, alpha = .4, size = 1) +
    stat_summary(aes(group = unknown_fraction), fun = median, geom = 'point', colour = '#CC6677', size = 3) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal() +
    xlab('True tumor fraction') +
    ylab('Unknown proportion estimated by SECRET') +
    labs(title = 'SECRET: estimated unknown content', subtitle = 'points: samples; red: median per tumor fraction')
  ggsave(file.path(output_dir, 'unknown_content_secret_unknown.pdf'), secret.plot, width = 6, height = 6)
  write_csv(secret_unknown %>% group_by(unknown_fraction) %>%
              dplyr::summarize(median_estimated_unknown = median(predicted_value), n = n()),
            file.path(output_dir, 'unknown_content_secret_unknown.csv'))
}

cat('figures written to', output_dir, '\n')
