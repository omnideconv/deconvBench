#!/usr/bin/Rscript

library(SingleCellExperiment)
library(parallel)

"Usage: 
  preprocessSingleCellNF.R <sc_name> <sc_path> <method> <subset_value> <replicate> <preProcess_dir> <baseDir>
  
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<method> deconvolution method, for which the preprocessing is performed. Important to select the correctly normalized sc dataset
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<replicate> value of replicate number
<preProcess_dir> directory where subsampling datasets are stored
<baseDir> nextflow base directory" -> doc


args <- docopt::docopt(doc)

sc_dataset <- args$sc_name
sc_path <- args$sc_path
method <- args$method
baseDir <- args$baseDir

source(paste0(baseDir, '/bin/utils.R'))
method_normalizations <- read.table(paste0(baseDir, '/optimal_normalizations.csv'), sep = ',', header = TRUE)
sc_norm <- method_normalizations[method_normalizations$method == method, 2]
print(paste0('Method: ', method, '; sc-norm: ', sc_norm))

subset_value <- as.numeric(args$subset_value)
replicate <- as.numeric(args$replicate)

output_dir <- paste0(args$preProcess_dir,'/',sc_dataset,'_',sc_norm,'_perc',subset_value,'_rep',replicate)

# check if pre-processed files already exist to save time
if(dir.exists(output_dir)){
  # check if all files are present
  if(all(c('batch.rds','celltype_annotations.rds') %in% list.files(output_dir))){
    cat('Preprocessing with given parameters has already been done and will be skipped.')
    quit(save='no')
  }
}

# read scRNA-seq count matrix
if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_counts.rds'))
} else {
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_norm_counts.rds'))
}
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_dataset, 'celltype_annotations.rds'))
sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))

n_cell_types <- length(unique(sc_celltype_annotations))
n_cells_per_ct <- table(sc_celltype_annotations) 


cat(paste0("Preprocessing sc_matrix file with subset_value=", subset_value, " and replicate=", replicate, "\n"))

# build seed for subsampling:
# the goal is to select the cells based on the subset size, replicate id and single-cell dataset
# other parameters (like sc_norm) should not influence cell sampling, so that we have the same
# cells selected even though the normlization changes
# this only works correctly to subset sizes > 0.1% (which should be fine)
seed_sc_name <- char2seed(sc_dataset)
seed <- as.integer((seed_sc_name + replicate + subset_value*1000) %% (2^31-1))
print(seed)
set.seed(seed)

total_cells <- length(sc_celltype_annotations)
min_n_cells <- 20 # fix the minimal number of cells per celltype. If less are available, all cells are selected

anno_tbl <- sc_celltype_annotations |>
  as.tibble() |>
  mutate(index = 1:length(sc_celltype_annotations)) |>
  group_by(value)

celltype_counts <- anno_tbl %>%
  dplyr::count(value, name = "n")

anno_subset <- anno_tbl %>%
  inner_join(celltype_counts, by = "value") %>%
  group_by(value) %>%
  group_modify(~ {
    n_total <- nrow(.x)
    target_n1 = round(n_total / total_cells * subset_value)
    target_n2 = max(target_n1, min_n_cells)
    target_n_final = min(target_n2, n_total)
    slice_sample(.x, n = target_n_final, replace = FALSE)
  }) %>%
  ungroup()


subset_matrix <- sc_matrix[,anno_subset$index]
subset_annot <- sc_celltype_annotations[anno_subset$index]
subset_batch <- sc_batch[anno_subset$index]


dir.create(output_dir, recursive = T, showWarnings = TRUE)
print(output_dir)

if(sc_norm == 'counts'){
  saveRDS(subset_matrix, paste0(output_dir,'/matrix_counts.rds'))
}else{
  saveRDS(subset_matrix, paste0(output_dir,'/matrix_norm_counts.rds'))
}
saveRDS(subset_annot, paste0(output_dir,'/celltype_annotations.rds'))
saveRDS(subset_batch, paste0(output_dir,'/batch.rds'))

cat(paste0("Preprocessed sc_matrix file stored in '", output_dir, "'\n"))
