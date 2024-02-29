#!/usr/bin/Rscript

print("Starting analysis script for deconvolution of simulated data ...")

library(docopt)
library(Biobase)
library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)

"Usage:
  nalysisNF_mixed_simulations.R <sc_name> <sc_path> <bulk_name> <bulk_path> <preprocess_dir> <replicates> <deconv_method> <results_dir> <ncores> <baseDir>
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<bulk_name> name of simulated bulk RNAseq dataset
<bulk_path> path to simulated bulk datasets
<preprocess_dir> preprocessing directory where pseudo-bulks are stored
<replicates> number of replicates
<deconv_method>  deconv method
<results_dir> results (base) directory
<ncores> number of cores to use for method (if available)
<baseDir> nextflow base directory" -> doc


args <- docopt::docopt(doc)

sc_dataset <- args$sc_name
sc_path <- args$sc_path
bulk_name <- args$bulk_name
bulk_path <- args$bulk_path
preprocess_dir <- args$preprocess_dir
method <- args$deconv_method
res_base_path <- args$results_dir
ncores <- as.numeric(args$ncores) # in case a method can use multiple cores
replicates <- as.numeric(args$replicates)
baseDir <- args$baseDir

source(paste0(baseDir, '/bin/utils.R'))
method_normalizations <- read.table(paste0(baseDir, '/optimal_normalizations.csv'), sep = ',', header = TRUE)
# Here we need to filter for those cell types that are in the simulated dataset. 
common.cells <- readRDS(file.path(preprocess_dir, 'cell_types.rds'))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_dataset, 'celltype_annotations.rds'))
position_vector <- sc_celltype_annotations %in% common.cells
sc_celltype_annotations <- sc_celltype_annotations[position_vector]

sc_norm <- method_normalizations[method_normalizations$method == method, 2]
bulk_norm <- method_normalizations[method_normalizations$method == method, 3]

if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_counts.rds'))[, position_vector] 
} else {
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_norm_counts.rds'))[, position_vector] 
}

sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))[position_vector]

# Subsetting the cells
subset_list <- subset_cells(sc_matrix, sc_celltype_annotations, sc_batch, 200, 22)

sc_matrix <- subset_list$data
sc_celltype_annotations <- subset_list$annotations
sc_batch <- subset_list$batch_id

res_path_normal <- paste0(res_base_path, '/', method, '_', bulk_name)

# Docker needs that
res_path_normal <- tolower(res_path_normal)

dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)

# We need a bulk for the creation of the signature matrix, because some methods filter for the common genes between sc and bulk. 
# In this case it is not relevant to compute the signature n times for each replicate

bulk_matrix <- readRDS(file.path(bulk_path,  paste0('replicate_1'), paste0('simulation_', bulk_norm, '.rds')))
bulk_matrix <- as.matrix(bulk_matrix)

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

for(r in 1:replicates){

    res_path <- file.path(res_path_normal, paste0(bulk_name, '_pseudobulk_', sc_dataset, '_signature'), paste0('replicate_', r))
    dir.create(res_path, recursive = TRUE, showWarnings = TRUE)
    bulk_matrix <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0('simulation_', bulk_norm, '.rds')))
    bulk_matrix <- as.matrix(bulk_matrix)

    # Deconvolution
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

    true_fractions <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0('simulation_facs.rds')))

    results_list = list(
        'deconvolution' = deconvolution, 
        'true_cell_fractions' = true_fractions
    )

    if(method =='autogenes' | method == 'scaden'){
    unlink(signature)
    }

saveRDS(results_list, file=paste0(res_path, "/deconvolution.rds")) 

}
