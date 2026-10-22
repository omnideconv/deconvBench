# Missing cell types (as Fig. S7 of the omnideconv benchmark) with the results of the rectangle revision:
# HaoCleaned-sampled as reference, hoek and finotello as bulks, eight omnideconv methods and rectangle.
# Based on supplement/fig_s7.R; the difference of the RMSE with and without the missing cell type is
# the same, now for all methods and one reference.
#
# Usage: Rscript fig_missing_cell_types.R [<results_dir>] [<output_dir>]

library(tidyverse)
library(ggpubr)
library(stringi)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else '/nfs/data/omnideconv_benchmarking_clean/benchmark_results/rectangle_revision/missing_cell_types'
output_dir <- if (length(args) >= 2) args[2] else '/nfs/proj/omnideconv_benchmarking/omnideconv/rectangle_revision/figures'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_names <- c(autogenes = 'AutoGeneS', bayesprism = 'BayesPrism', bisque = 'Bisque', cibersortx = 'CIBERSORTx',
                  dwls = 'DWLS', music = 'MuSiC', rectangle = 'Rectangle', scaden = 'Scaden', scdc = 'SCDC')

# 1: List directories, methods, cell types

files <- list.files(results_dir, full.names = FALSE, recursive = TRUE)
files <- files[grepl('deconvolution.rds', files, fixed = TRUE)]

metadata.table <- tibble(path = files,
                         method = map_chr(files, function(x) strsplit(x, split = '/')[[1]][1]),
                         metadata = map_chr(files, function(x) strsplit(x, split = '/')[[1]][2])) %>%
  mutate(metadata = sub('_no_', '_', metadata)) %>%
  mutate(sc_dataset = map_chr(metadata, function(x) stri_split_fixed(x, pattern = '_', n = 3)[[1]][1]),
         bulk_dataset = map_chr(metadata, function(x) stri_split_fixed(x, pattern = '_', n = 3)[[1]][2]),
         missing_celltype = gsub('_', ' ', map_chr(metadata, function(x) stri_split_fixed(x, pattern = '_', n = 3)[[1]][3]))) %>%
  select(-metadata)

cat('result files per method:\n'); print(table(metadata.table$method))

# 2: Correlation and RMSE per predicted cell type, with and without each missing cell type

read_metrics <- function(i) {
  res <- readRDS(file.path(results_dir, metadata.table$path[i]))
  corr <- as.data.frame(res$cor_cell_type) %>% select(cell_type, cor)
  rmse <- as.data.frame(res$rmse_cell_type) %>% select(cell_type, RMSE)
  left_join(corr, rmse, by = 'cell_type') %>%
    mutate(cell_type = gsub('xxxx', ' ', cell_type)) %>%
    bind_cols(metadata.table[i, c('method', 'sc_dataset', 'bulk_dataset', 'missing_celltype')])
}
metrics.data <- bind_rows(lapply(seq_len(nrow(metadata.table)), read_metrics))
metrics.data$missing_celltype[metrics.data$missing_celltype == 'all-cells'] <- 'None'

celltypes_ordered_similarity <- c('Tregs', 'T cells CD4 conv', 'T cells CD8', 'ILC', 'Plasma cells', 'NK cells',
                                  'B cells', 'Monocytes', 'mDC', 'pDC', 'Platelet')

metrics.all.celltypes <- metrics.data %>%
  filter(missing_celltype == 'None') %>%
  select(cell_type, method, sc_dataset, bulk_dataset, RMSE_all_celltypes = RMSE)
metrics.plot <- metrics.data %>%
  filter(missing_celltype != 'None') %>%
  select(cell_type, method, sc_dataset, bulk_dataset, missing_celltype, RMSE_defective = RMSE) %>%
  left_join(metrics.all.celltypes, by = c('cell_type', 'method', 'sc_dataset', 'bulk_dataset')) %>%
  mutate(abs_diff_rmse = RMSE_defective - RMSE_all_celltypes,
         method_name = factor(method_names[method], levels = method_names))

write_csv(metrics.plot, file.path(output_dir, 'missing_cell_types_rmse_difference.csv'))

# 3: Heatmaps: RMSE with the missing cell type minus RMSE with the complete reference

heatmap_plot <- function(bulk, title) {
  ggplot(metrics.plot[metrics.plot$bulk_dataset == bulk, ], aes(x = missing_celltype, y = cell_type, fill = abs_diff_rmse)) +
    geom_tile() +
    coord_equal() +
    geom_text(aes(label = round(abs_diff_rmse, 2)), size = 2) +
    facet_wrap(. ~ method_name, ncol = 3) +
    theme_minimal() +
    scale_fill_gradient2(mid = 'white', high = '#8a7495', low = '#7f9574', name = 'RMSE defective - \nRMSE full') +
    rotate_x_text(angle = 60) +
    scale_x_discrete(limits = celltypes_ordered_similarity) +
    labs(y = 'Predicted cell type\n', x = 'Missing cell type', title = title)
}

for (ext in c('pdf', 'png')) {
  ggsave(file.path(output_dir, paste0('missing_cell_types_finotello.', ext)), heatmap_plot('finotello', 'Finotello dataset'), dpi = 200, width = 12, height = 12, bg = 'white')
  ggsave(file.path(output_dir, paste0('missing_cell_types_hoek.', ext)), heatmap_plot('hoek', 'Hoek dataset'), dpi = 200, width = 12, height = 12, bg = 'white')
}

cat('figures written to', output_dir, '\n')
