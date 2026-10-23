# Resolution analysis of the rectangle revision: bump chart of the method ranks by RMSE for every cell
# type of the fine annotation (rank 1 = lowest RMSE), one chart per dataset (wu, lambrechts).
# The RMSE of a cell type is computed over all pseudo-bulk samples of the five replicates.
#
# Usage: Rscript fig_resolution_bump.R [<results_dir>] [<output_dir>]

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
results_dir <- if (length(args) >= 1) args[1] else '/nfs/data/omnideconv_benchmarking_clean/benchmark_results/rectangle_revision/resolution'
output_dir <- if (length(args) >= 2) args[2] else '/nfs/proj/omnideconv_benchmarking/omnideconv/rectangle_revision/figures'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_names <- c(rectangle = 'Rectangle', scdc = 'SCDC', autogenes = 'AutoGeneS', music = 'MuSiC', `music-hc` = 'MuSiC-hc',
                  scaden = 'Scaden', bayesprism = 'BayesPrism', bisque = 'Bisque', cibersortx = 'CIBERSORTx', dwls = 'DWLS',
                  secret = 'SECRET')
method_colors <- c(Rectangle = '#000000', SCDC = '#AFAFAF', AutoGeneS = '#7BC8A4', MuSiC = '#F2D649', `MuSiC-hc` = '#A39B2F',
                   Scaden = '#DDBF94', BayesPrism = '#EE8C5E', Bisque = '#8E9CC8', CIBERSORTx = '#DE89C1', DWLS = '#4F8FC0',
                   SECRET = '#B0413E')

# The methods write the cell-type names in different forms: rectangle with underscores instead of
# spaces, Scaden with dots for every other character, CIBERSORTx with 'xxxx' for spaces. They are
# matched to the names of the ground truth on letters and digits only.
name_key <- function(x) tolower(gsub('[^A-Za-z0-9]', '', gsub('xxxx', '', x)))

read_dataset <- function(dataset) {
  rows <- list()
  for (method_dir in list.dirs(results_dir, recursive = FALSE)) {
    method <- basename(method_dir)
    files <- list.files(file.path(method_dir, paste0(dataset, '_resolution_analysis_sim_fine_annot')),
                        pattern = 'deconvolution.rds', recursive = TRUE, full.names = TRUE)
    for (f in files) {
      res <- readRDS(f)
      estimated <- as.matrix(res$deconvolution)
      truth <- as.matrix(res$true_cell_fractions)              # cell types x samples
      keys <- name_key(rownames(truth))
      stopifnot(!any(duplicated(keys)))
      matched <- match(keys, name_key(colnames(estimated)))
      # a cell type without an estimate (e.g. left out by MuSiC for lack of donors) counts as 0
      est <- sapply(seq_along(keys), function(k) if (is.na(matched[k])) rep(0, nrow(estimated)) else estimated[, matched[k]])
      colnames(est) <- rownames(truth)
      rownames(est) <- rownames(estimated)
      rows[[length(rows) + 1]] <- tibble(
        method = method, replicate = basename(dirname(f)),
        sample = rep(rownames(est), times = ncol(est)), cell_type = rep(colnames(est), each = nrow(est)),
        estimated = as.vector(est), true = as.vector(t(truth[, rownames(est)])),
        has_estimate = rep(!is.na(matched), each = nrow(est)))
    }
  }
  bind_rows(rows)
}

bump_chart <- function(dataset) {
  data <- read_dataset(dataset)
  coverage <- data %>% group_by(method) %>% summarise(cell_types = n_distinct(cell_type),
                                                      without_estimate = n_distinct(cell_type[!has_estimate]),
                                                      replicates = n_distinct(replicate), .groups = 'drop')
  cat('\n', dataset, ':\n', sep = ''); print(coverage)

  ranks <- data %>%
    group_by(method, cell_type) %>%
    summarise(rmse = sqrt(mean((estimated - true)^2)), .groups = 'drop') %>%
    group_by(cell_type) %>%
    mutate(rank = rank(rmse, ties.method = 'min')) %>%
    ungroup() %>%
    mutate(method_name = factor(method_names[method], levels = method_names))
  write_csv(ranks, file.path(output_dir, paste0('resolution_', dataset, '_fine_rmse_ranks.csv')))

  n_methods <- n_distinct(ranks$method)
  n_types <- n_distinct(ranks$cell_type)
  ranks <- ranks %>% mutate(is_rectangle = method == 'rectangle')
  plot <- ggplot(ranks, aes(x = cell_type, y = rank, group = method_name, color = method_name)) +
    geom_line(data = filter(ranks, !is_rectangle), linewidth = 0.7, alpha = 0.85) +
    geom_point(data = filter(ranks, !is_rectangle), size = 2) +
    geom_line(data = filter(ranks, is_rectangle), linewidth = 1.3) +
    geom_point(data = filter(ranks, is_rectangle), size = 2.8) +
    scale_y_reverse(breaks = seq_len(n_methods), expand = expansion(add = 0.4)) +
    scale_color_manual(values = method_colors, drop = TRUE) +
    guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2, size = 2.5))) +
    labs(x = NULL, y = 'RMSE rank (1 = best)', color = NULL,
         title = paste0(tools::toTitleCase(dataset), ': method ranks by RMSE across ', n_types, ' fine-grained cell types')) +
    theme_classic(base_size = 11) +
    theme(legend.position = 'top',
          panel.grid.major.y = element_line(colour = 'grey90'),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  width <- max(8, 0.3 * n_types + 2)
  ggsave(file.path(output_dir, paste0('resolution_', dataset, '_fine_bump.pdf')), plot, width = width, height = 4.8)
  ggsave(file.path(output_dir, paste0('resolution_', dataset, '_fine_bump.png')), plot, width = width, height = 4.8, dpi = 200)

  summary <- ranks %>% group_by(method_name) %>%
    summarise(mean_rank = mean(rank), times_rank_1 = sum(rank == 1), median_rmse = median(rmse), .groups = 'drop') %>%
    arrange(mean_rank)
  write_csv(summary, file.path(output_dir, paste0('resolution_', dataset, '_fine_rank_summary.csv')))
  print(summary)
}

for (dataset in c('wu', 'lambrechts')) bump_chart(dataset)
cat('figures written to', output_dir, '\n')
