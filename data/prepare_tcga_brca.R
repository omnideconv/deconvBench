# Input of the TCGA-BRCA deconvolution: the bulk dataset tcga-brca and the reference wu-normal-tcga.
#
# Reference: wu at the normal annotation level without cycling, myoepithelial, luminal and cancer
# epithelial cell types, and without the genes that match ^RP[SL]|^RP\d|^RPSA. The CPM matrix is
# computed again from the remaining counts. With <max_cells> > 0, at most that many cells per cell
# type are kept (subset_cells of the pipeline, seed 22, as in the analysis workflows) and the
# reference is wu-normal-tcga; with 0 all cells are kept and the reference is wu-normal-tcga-full.
# Bulk: the TCGA counts and TPM as they are, as gene x sample matrices.
#
# Usage: Rscript prepare_tcga_brca.R <max_cells> [<data_dir>] [<pipeline_dir>]

suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
max_cells <- as.numeric(args[1])
data_dir <- if (length(args) >= 2) args[2] else '/nfs/data/omnideconv_benchmarking_clean/data'
pipeline_dir <- if (length(args) >= 3) args[3] else '/nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/pipeline'

excluded_cell_types <- c('Cycling T-cells', 'Cycling_Myeloid', 'Cancer Cycling', 'Cycling PVL', 'Myoepithelial',
                         'Luminal Progenitors', 'Mature Luminal', 'Cancer Her2 SC', 'Cancer LumB SC',
                         'Cancer Basal SC', 'Cancer LumA SC')
excluded_genes <- '^RP[SL]|^RP\\d|^RPSA'

#### reference ####
wu <- file.path(data_dir, 'singleCell', 'wu')
counts <- readRDS(file.path(wu, 'matrix_counts.rds'))
annotations <- as.character(readRDS(file.path(wu, 'celltype_annotations_normal.rds')))
batch <- as.character(readRDS(file.path(wu, 'batch.rds')))
stopifnot(all(excluded_cell_types %in% annotations), ncol(counts) == length(annotations))

keep_cells <- !annotations %in% excluded_cell_types
keep_genes <- !grepl(excluded_genes, rownames(counts))
counts <- counts[keep_genes, keep_cells]
annotations <- annotations[keep_cells]
batch <- batch[keep_cells]
cat('cells:', ncol(counts), '; cell types:', length(unique(annotations)), '; genes:', nrow(counts),
    '(', sum(!keep_genes), 'removed )\n')

if (max_cells > 0) {
  source(file.path(pipeline_dir, 'bin', 'utils.R'))
  subset <- subset_cells(counts, annotations, batch, max_cells, 22)
  counts <- as.matrix(subset$data)
  annotations <- as.character(subset$annotations)
  batch <- as.character(subset$batch_id)
  cat('after downsampling to at most', max_cells, 'cells per cell type:', ncol(counts), 'cells\n')
}
print(table(annotations))

cpm <- t(t(counts) / colSums(counts)) * 1e6

out <- file.path(data_dir, 'singleCell', if (max_cells > 0) 'wu-normal-tcga' else 'wu-normal-tcga-full')
dir.create(out, showWarnings = FALSE)
saveRDS(counts, file.path(out, 'matrix_counts.rds'))
saveRDS(cpm, file.path(out, 'matrix_norm_counts.rds'))
saveRDS(annotations, file.path(out, 'celltype_annotations.rds'))
saveRDS(batch, file.path(out, 'batch.rds'))
writeLines(c('wu at the normal annotation level, prepared by benchmark/data/prepare_tcga_brca.R',
             paste('excluded cell types:', paste(excluded_cell_types, collapse = ', ')),
             paste('excluded genes (regex):', excluded_genes),
             paste('at most cells per cell type (0 = all):', max_cells, '(seed 22)'),
             'matrix_norm_counts.rds: CPM of the remaining genes'),
           file.path(out, 'README.txt'))

#### bulk ####
read_tcga <- function(f) as.matrix(read.table(file.path(data_dir, 'TCGA-BRCA', f), header = TRUE, sep = '\t', check.names = FALSE))
bulk_counts <- read_tcga('TCGA_counts_BRCA.txt')
bulk_tpm <- read_tcga('TCGA_tpm_BRCA.txt')
stopifnot(identical(dimnames(bulk_counts), dimnames(bulk_tpm)), !anyNA(bulk_counts), !anyNA(bulk_tpm))
bulk_out <- file.path(data_dir, 'TCGA-BRCA', 'tcga-brca')
dir.create(bulk_out, showWarnings = FALSE)
saveRDS(bulk_counts, file.path(bulk_out, 'tcga-brca_counts.rds'))
saveRDS(bulk_tpm, file.path(bulk_out, 'tcga-brca_tpm.rds'))
cat('bulk:', nrow(bulk_tpm), 'genes x', ncol(bulk_tpm), 'samples; genes shared with the reference:',
    length(intersect(rownames(bulk_tpm), rownames(counts))), '\n')
