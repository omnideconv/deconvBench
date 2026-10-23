# Description of the single-cell references of the rectangle revision: cells, genes and donors,
# cell types per annotation level with their hierarchy, rare cell types, detected genes and counts
# per cell type, and the cell-type composition of every donor.
#
# Usage: Rscript fig_reference_qc.R [<output_dir>]

suppressMessages({library(tidyverse); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[1] else '/nfs/proj/omnideconv_benchmarking/omnideconv/rectangle_revision/figures/reference'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- '/nfs/data/omnideconv_benchmarking_clean/data'
datasets <- c(lambrechts = 'singleCell/lambrechts', wu = 'singleCell/wu', `HaoCleaned-sampled` = 'singleCellXGeneCensus/HaoCleaned-sampled')
rare_threshold <- 50

summary_rows <- list(); level_rows <- list(); hierarchy_rows <- list(); qc_rows <- list(); mapping_rows <- list()

for (ds in names(datasets)) {
  dir <- file.path(data_dir, datasets[[ds]])
  counts <- readRDS(file.path(dir, 'matrix_counts.rds'))
  donor <- as.character(readRDS(file.path(dir, if (file.exists(file.path(dir, 'donor.rds'))) 'donor.rds' else 'batch.rds')))
  levels <- c(coarse = 'celltype_annotations_coarse.rds', normal = 'celltype_annotations_normal.rds', fine = 'celltype_annotations_fine.rds')
  levels <- levels[file.exists(file.path(dir, levels))]
  if (length(levels) == 0) levels <- c(normal = 'celltype_annotations.rds')
  ann <- lapply(levels, function(f) as.character(readRDS(file.path(dir, f))))
  stopifnot(all(lengths(ann) == ncol(counts)), length(donor) == ncol(counts))

  # detected genes and counts per cell, in blocks to keep the dense matrices manageable
  blocks <- split(seq_len(ncol(counts)), ceiling(seq_len(ncol(counts)) / 5000))
  n_genes <- unlist(lapply(blocks, function(b) colSums(as.matrix(counts[, b]) > 0)))
  n_counts <- unlist(lapply(blocks, function(b) colSums(as.matrix(counts[, b]))))

  finest <- names(ann)[length(ann)]
  coarsest <- names(ann)[1]
  qc_rows[[ds]] <- tibble(dataset = ds, cell_type = ann[[finest]], group = ann[[coarsest]], donor = donor,
                          n_genes = n_genes, n_counts = n_counts)

  summary_rows[[ds]] <- tibble(dataset = ds, cells = ncol(counts), genes = nrow(counts), donors = length(unique(donor)),
                               cells_per_donor_min = min(table(donor)), cells_per_donor_median = median(table(donor)),
                               cells_per_donor_max = max(table(donor)),
                               levels = paste(names(ann), sapply(ann, function(a) length(unique(a))), sep = ': ', collapse = ', '),
                               median_genes_per_cell = median(n_genes), median_counts_per_cell = median(n_counts))

  for (lvl in names(ann)) {
    level_rows[[paste(ds, lvl)]] <- tibble(dataset = ds, level = lvl, cell_type = ann[[lvl]], donor = donor) %>%
      group_by(dataset, level, cell_type) %>%
      summarise(cells = n(), donors = n_distinct(donor), max_donor_share = max(table(donor)) / n(), .groups = 'drop')
  }
  hierarchy_rows[[ds]] <- as_tibble(ann) %>% mutate(donor = donor) %>%
    group_by(across(all_of(names(ann)))) %>% summarise(cells = n(), donors = n_distinct(donor), .groups = 'drop') %>%
    mutate(dataset = ds, .before = 1)

  # the mapping file that the resolution simulation uses, against the annotations of the cells
  mapping_file <- file.path(dir, 'cell_type_mappings.csv')
  if (file.exists(mapping_file) && 'fine' %in% names(ann)) {
    m <- read.csv(mapping_file, check.names = FALSE)
    observed <- unique(tibble(Fine = ann$fine, Normal = ann$normal, Coarse = ann$coarse))
    joined <- full_join(observed, m, by = 'Fine', suffix = c('', '_mapping'))
    mapping_rows[[ds]] <- joined %>%
      filter(is.na(Normal) | is.na(Normal_mapping) | Normal != Normal_mapping | Coarse != Coarse_mapping) %>%
      mutate(dataset = ds, .before = 1)
  }
  rm(counts); gc()
}

summary_table <- bind_rows(summary_rows)
level_table <- bind_rows(level_rows)
hierarchy <- bind_rows(hierarchy_rows)
qc <- bind_rows(qc_rows)
mapping_issues <- bind_rows(mapping_rows)

write_csv(summary_table, file.path(output_dir, 'reference_summary.csv'))
write_csv(level_table, file.path(output_dir, 'reference_cell_types.csv'))
write_csv(hierarchy, file.path(output_dir, 'reference_hierarchy.csv'))
write_csv(mapping_issues, file.path(output_dir, 'reference_mapping_issues.csv'))
print(summary_table, width = 200)
cat('\nrare cell types (<', rare_threshold, 'cells):\n')
print(level_table %>% filter(cells < rare_threshold) %>% arrange(dataset, level, cells), n = 100)
cat('\nmapping file vs annotations:\n'); print(mapping_issues, width = 200)

# detected genes per cell, per finest cell type, faceted by the coarsest group
for (ds in names(datasets)) {
  d <- qc %>% filter(dataset == ds)
  type_order <- d %>% group_by(group, cell_type) %>% summarise(m = median(n_genes), .groups = 'drop') %>% arrange(group, m)
  d <- d %>% mutate(cell_type = factor(cell_type, levels = unique(type_order$cell_type)))
  n_types <- n_distinct(d$cell_type)
  one_level <- n_distinct(d$group) == n_types
  p <- ggplot(d, aes(x = n_genes, y = cell_type)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3, fill = '#C9D9DC', colour = '#2F5E66', linewidth = 0.3, width = 0.7) +
    scale_x_log10() +
    labs(x = 'Detected genes per cell (log scale)', y = NULL,
         title = paste0(ds, ': detected genes per cell and cell type'),
         subtitle = paste0(nrow(d), ' cells; boxes: median and quartiles')) +
    theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0, hjust = 0))
  if (!one_level) p <- p + facet_grid(group ~ ., scales = 'free_y', space = 'free_y')
  ggsave(file.path(output_dir, paste0('qc_genes_', ds, '.png')), p, width = 8, height = 1.2 + 0.17 * n_types + (if (one_level) 0 else 1), dpi = 200, bg = 'white')
  ggsave(file.path(output_dir, paste0('qc_genes_', ds, '.pdf')), p, width = 8, height = 1.2 + 0.17 * n_types + (if (one_level) 0 else 1))
}

# cell-type composition of every donor at the coarsest level
comp <- qc %>% count(dataset, donor, group) %>% group_by(dataset, donor) %>% mutate(share = n / sum(n), cells = sum(n)) %>% ungroup()
for (ds in names(datasets)) {
  d <- comp %>% filter(dataset == ds)
  groups <- sort(unique(d$group))
  palette <- setNames(colorRampPalette(c('#1F6F78', '#8FB8A8', '#E0C878', '#C7704F', '#7A5C8E', '#4F6D9A'))(length(groups)), groups)
  donor_order <- d %>% group_by(donor) %>% summarise(cells = first(cells)) %>% arrange(cells) %>% pull(donor)
  p <- ggplot(d, aes(x = share, y = factor(donor, levels = donor_order), fill = group)) +
    geom_col(width = 0.8, colour = 'white', linewidth = 0.2) +
    geom_text(data = distinct(d, donor, cells), aes(x = 1.01, y = donor, label = format(cells, big.mark = ',')),
              inherit.aes = FALSE, hjust = 0, size = 2.6, colour = 'grey30') +
    scale_x_continuous(labels = scales::percent, limits = c(0, 1.12), breaks = seq(0, 1, 0.25)) +
    scale_fill_manual(values = palette) +
    labs(x = 'Share of the donor\'s cells', y = NULL, fill = NULL, title = paste0(ds, ': cell types per donor'),
         subtitle = 'Coarsest annotation level; number of cells on the right') +
    theme_minimal(base_size = 10) + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = 'bottom') +
    guides(fill = guide_legend(nrow = 2))
  h <- 1.8 + 0.2 * length(donor_order)
  ggsave(file.path(output_dir, paste0('donor_composition_', ds, '.png')), p, width = 8, height = h, dpi = 200, bg = 'white')
  ggsave(file.path(output_dir, paste0('donor_composition_', ds, '.pdf')), p, width = 8, height = h)
}
cat('written to', output_dir, '\n')
