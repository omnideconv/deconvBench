#!/usr/bin/env Rscript
# Build the CI test data in data/ci/ from the small datasets in data/. Runs inside the standard
# image (needs SimBu). The data contain the edge cases that broke runs on real data:
#   - annotations stored as a factor with an unused level (HaoCleaned broke AutoGeneS and SCDC)
#   - a cell type whose cells all come from one donor (MuSiC leaves it out, music-hc failed on it;
#     a split by donor cannot put it into both parts)
#   - cell-type names with '+', '-', '.', '_' and spaces (Scaden, rectangle and CIBERSORTx change them)
#   - two cell types whose names contain each other (broke the reference of SECRET)
#   - a bulk dataset without ground truth and with non-integer counts (TCGA)
#   - a bulk named hoek (its ground truth sums T cell subtypes)
#
# Usage: Rscript make_ci_fixture.R <repo_data_dir> <out_dir>

suppressMessages({library(SimBu); library(Matrix)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: make_ci_fixture.R <repo_data_dir> <out_dir>')
data_dir <- args[1]
out_dir <- args[2]
set.seed(1)

sc_in <- file.path(data_dir, 'single_cell', 'sc1')
counts <- readRDS(file.path(sc_in, 'matrix_counts.rds'))
annot <- as.character(readRDS(file.path(sc_in, 'celltype_annotations.rds')))
coarse <- as.character(readRDS(file.path(sc_in, 'celltype_annotations_coarse.rds')))
batch <- as.character(readRDS(file.path(sc_in, 'batch.rds')))

# cell types: PBMC names kept where ground truths and the hoek sum need them
keep <- annot != 'ILC'
counts <- counts[, keep]; annot <- annot[keep]; coarse <- coarse[keep]; batch <- batch[keep]
mono <- which(annot == 'Monocytes')
annot[mono[seq(2, length(mono), 2)]] <- 'Monocytes classical'
annot[annot == 'pDC'] <- 'pDC_c4-IRF7+'
batch[annot == 'pDC_c4-IRF7+'] <- 'P7'
annot[annot == 'Platelet'] <- 'Platelets.CD41+'

fine <- annot
for (ct in c('B cells', 'T cells CD8')) {
  pos <- which(annot == ct)
  fine[pos] <- paste0(ct, c(' s1', ' s2'))[(seq_along(pos) %% 2) + 1]
}

# 3,000 genes with the highest variance of log CPM
cpm_all <- t(t(as.matrix(counts)) / colSums(counts)) * 1e6
genes <- names(sort(apply(log1p(cpm_all), 1, var), decreasing = TRUE))[1:3000]
counts <- as.matrix(counts[genes, ])
cpm <- t(t(counts) / colSums(counts)) * 1e6

write_sc <- function(name, idx, with_coarse = TRUE, with_resolution = TRUE) {
  out <- file.path(out_dir, 'single_cell', name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(counts[, idx], file.path(out, 'matrix_counts.rds'))
  saveRDS(cpm[, idx], file.path(out, 'matrix_norm_counts.rds'))
  # a factor with an unused level, as in HaoCleaned
  saveRDS(factor(annot[idx], levels = c(sort(unique(annot[idx])), 'unused type')), file.path(out, 'celltype_annotations.rds'))
  saveRDS(batch[idx], file.path(out, 'batch.rds'))
  simbu_ds <- SimBu::dataset(annotation = data.frame(ID = colnames(counts)[idx], cell_type = annot[idx]),
                             count_matrix = Matrix(counts[, idx], sparse = TRUE),
                             tpm_matrix = Matrix(cpm[, idx], sparse = TRUE), name = name)
  saveRDS(simbu_ds, file.path(out, 'simbu_ds.rds'))
  if (with_resolution) {
    saveRDS(annot[idx], file.path(out, 'celltype_annotations_normal.rds'))
    saveRDS(fine[idx], file.path(out, 'celltype_annotations_fine.rds'))
    mapping <- unique(data.frame(Fine = fine[idx], Normal = annot[idx], Coarse = coarse[idx]))
    write.csv(mapping[order(mapping$Fine), ], file.path(out, 'cell_type_mappings.csv'), row.names = FALSE, quote = FALSE)
    if (with_coarse) saveRDS(coarse[idx], file.path(out, 'celltype_annotations_coarse.rds'))
  }
}
all_cells <- seq_along(annot)
write_sc('sc-ci', all_cells)
write_sc('sc-ci-nc', all_cells, with_coarse = FALSE)
write_sc('sc-ci2', which(all_cells %% 3 != 0 & !annot %in% c('Platelets.CD41+', 'Plasma cells')), with_resolution = FALSE)

# bulks
bulk_dir <- function(name) { d <- file.path(out_dir, 'bulks', name); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
for (name in c('bulk-ci', 'hoek')) {
  d <- bulk_dir(name)
  for (suffix in c('counts', 'tpm', 'facs')) {
    saveRDS(readRDS(file.path(data_dir, 'bulks', 'bulk1', paste0('bulk1_', suffix, '.rds'))), file.path(d, paste0(name, '_', suffix, '.rds')))
  }
}
d <- bulk_dir('bulk-nofacs')
saveRDS(as.matrix(readRDS(file.path(data_dir, 'bulks', 'bulk2', 'bulk2_counts.rds'))) * 1.37, file.path(d, 'bulk-nofacs_counts.rds'))
saveRDS(readRDS(file.path(data_dir, 'bulks', 'bulk2', 'bulk2_tpm.rds')), file.path(d, 'bulk-nofacs_tpm.rds'))

cat('CI data written to', out_dir, '\n')
print(table(annot, batch))
