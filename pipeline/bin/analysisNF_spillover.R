#!/usr/bin/Rscript

print("Starting analysis script [spillover analysis] ...")

library(omnideconv)
reticulate::use_condaenv(condaenv = "r-omnideconv", required = TRUE)
sessionInfo()
reticulate::py_config()

"Usage:
  analysisNF_spillover.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <cell_types> <results_dir> <ncores> <baseDir>
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<bulk_name> name of simulated bulk RNAseq dataset
<bulk_path> path to simulated bulk datasets
<deconv_method>  deconv method
<cell_types> cell types used in the simulation
<results_dir> results (base) directory
<ncores> number of cores to use for method (if available)
<baseDir> nextflow base directory" -> doc


args <- docopt::docopt(doc)

sc_dataset <- args$sc_name
sc_path <- args$sc_path
bulk_name <- args$bulk_name
bulk_path <- args$bulk_path
method <- args$deconv_method
res_base_path <- args$results_dir
ncores <- as.numeric(args$ncores) # in case a method can use multiple cores
baseDir <- args$baseDir

source(paste0(baseDir, '/bin/utils.R'))
method_normalizations <- read.table(paste0(baseDir, '/optimal_normalizations.csv'), sep = ',', header = TRUE)
sc_norm <- method_normalizations[method_normalizations$method == method, 2]
bulk_norm <- method_normalizations[method_normalizations$method == method, 3]
print(paste0('Method: ', method, '; sc-norm: ', sc_norm, '; bulk-norm: ', bulk_norm))

if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_counts.rds'))
} else {
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_norm_counts.rds'))
}

sc_celltype_annotations <- readRDS(file.path(sc_path, sc_dataset, 'celltype_annotations.rds'))
sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))

# Here we need to filter for those cell types that are in the simulated dataset. 
cat(paste0('Number of cells in sc dataset before filtering for target types: ', length(sc_celltype_annotations), '\n'))

cell_types_simulation <- gsub('\\[|]', '', args$cell_types)
cell_types_simulation <- strsplit(cell_types_simulation, ", ")[[1]]

position_vector <- (sc_celltype_annotations %in% cell_types_simulation)
sc_matrix <- sc_matrix[, position_vector]
sc_batch <- sc_batch[position_vector]
sc_celltype_annotations <- sc_celltype_annotations[position_vector]

cat(paste0('Number of cells in sc dataset before filtering for target types: ', length(sc_celltype_annotations), '\n'))

bulk_matrix <- as.matrix(readRDS(file.path(bulk_path, paste0(bulk_name, '_', bulk_norm, '.rds'))))

res_path_normal <- paste0(res_base_path, '/', method, "_", sc_dataset, '_', bulk_name)

dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)


signature <- signature_workflow_general(
  sc_matrix, 
  sc_celltype_annotations, 
  'normal', 
  sc_dataset, 
  sc_norm, 
  sc_batch, 
  method, 
  bulk_matrix,
  bulk_name, 
  bulk_norm, 
  ncores, 
  res_path_normal,
  baseDir=baseDir
)


deconvolution <- deconvolution_workflow_general(
  sc_matrix, 
  sc_celltype_annotations,
  'normal', 
  sc_dataset, 
  sc_norm, 
  sc_batch, 
  signature, 
  method, 
  bulk_matrix, 
  bulk_name, 
  bulk_norm, 
  ncores, 
  res_path_normal,
  baseDir=baseDir
)

true_fractions <- readRDS(file.path(bulk_path, paste0(bulk_name, '_facs.rds')))

results_list = list(
  'deconvolution' = deconvolution, 
  'true_cell_fractions' = true_fractions
)

saveRDS(results_list, file=paste0(res_path_normal, "/deconvolution_spillover.rds"))