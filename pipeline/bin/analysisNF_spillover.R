#!/usr/bin/Rscript

print("Starting analysis script [spillover analysis] ...")

library(docopt)
library(Biobase)
library(omnideconv)
library(tidyverse)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)

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

if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_counts.rds'))
} else {
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_norm_counts.rds'))
}

sc_celltype_annotations <- readRDS(file.path(sc_path, sc_dataset, 'celltype_annotations.rds'))
sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))


# Here we need to filter for those cell types that are in the simulated dataset. 

cell_types_simulation <- gsub('\\[|]', '', args$cell_types)
cell_types_simulation <- strsplit(cell_types_simulation, ",")[[1]]
print(cell_types_simulation)

position_vector <- sc_celltype_annotations %in% cell_types_simulation
sc_matrix <- sc_matrix[, position_vector]
sc_batch <- sc_batch[position_vector]
sc_celltype_annotations <- sc_celltype_annotations[position_vector]

bulk_matrix <- readRDS(file.path(bulk_path, paste0(bulk_name, '_', bulk_norm, '.rds')))
bulk_matrix <- as.matrix(bulk_matrix)

res_path_normal <- paste0(res_base_path, '/', method, "_", sc_dataset, '_', bulk_name)

dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)

subset_list <- subset_cells(sc_matrix, sc_celltype_annotations, sc_batch, 500, 22)

sc_matrix <- subset_list$data
sc_celltype_annotations <- subset_list$annotations
sc_batch <- subset_list$batch_id

# Signature building can be done here already, as the bulk dataset has no influence on signatures
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
  res_path_normal
)

# Deconvolution
for(cur_cell_type in cell_types_simulation){
  
  cur_cell_type <- gsub(' ', '_', cur_cell_type)
  bulk_matrix <- readRDS(file.path(bulk_path, paste0(bulk_name, '_', cur_cell_type, '_', bulk_norm, '.rds')))
  bulk_matrix <- as.matrix(bulk_matrix)
  
  true_fractions <- readRDS(file.path(bulk_path, paste0(bulk_name, '_', cur_cell_type, '_facs.rds')))

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
    res_path_normal
  )
  
  results_list = list(
    'deconvolution' = deconvolution, 
    'true_cell_fractions' = true_fractions
  )
  saveRDS(results_list, file=paste0(res_path_normal, "/deconvolution_spillover_", cur_cell_type, ".rds"))
  
}

if(method =='autogenes' | method == 'scaden'){
  unlink(signature)
}
