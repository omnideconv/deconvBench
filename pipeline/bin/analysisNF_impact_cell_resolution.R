#!/usr/bin/Rscript

print("Starting analysis script [impact cell resolution] ...")

library(docopt)
library(Biobase)
library(omnideconv)
library(tidyverse)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)

"Usage:
  analysisNF_impact_cell_resolution.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <cell_types_fine> <replicates> <results_dir> <ncores> <baseDir>
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<bulk_name> name of simulated bulk RNAseq dataset
<bulk_path> path to simulated bulk datasets
<deconv_method>  deconv method
<cell_types_fine> cell types used in the simulation
<replicates> number of replicates
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
replicates <- as.numeric(args$replicates)
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

annotation_path <- file.path(sc_path, sc_dataset, 'celltype_annotations_normal.rds')
sc_celltype_annotations <- readRDS(annotation_path)
sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))

sc_celltype_annotations_fine <- readRDS(file.path(gsub('celltype_annotations_normal.rds', 'celltype_annotations_fine.rds', annotation_path)))
sc_celltype_annotations_coarse <- readRDS(file.path(gsub('celltype_annotations_normal.rds', 'celltype_annotations_coarse.rds', annotation_path)))

# Here we need to filter for those cell types that are in the simulated dataset. 
# NOTE: in the resolution analysis the cell types are specified in terms of finest cell types

cell_types_fine <- gsub('\\[|]', '', args$cell_types_fine)
cell_types_fine <- strsplit(cell_types_fine, ",")[[1]]
print(cell_types_fine)

position_vector <- sc_celltype_annotations_fine %in% cell_types_fine
sc_matrix <- sc_matrix[, position_vector]
sc_batch <- sc_batch[position_vector]
sc_celltype_annotations_fine <- sc_celltype_annotations_fine[position_vector]
sc_celltype_annotations_coarse <- sc_celltype_annotations_coarse[position_vector]
sc_celltype_annotations <- sc_celltype_annotations[position_vector]

bulk_matrix <- readRDS(file.path(bulk_path, paste0('replicate_1'), paste0(bulk_name, '_', bulk_norm, '.rds'))) %>% as.matrix(.)

res_path_normal <- paste0(res_base_path, '/', method, '/', bulk_name, "_normal_annot")
res_path_fine <- paste0(res_base_path, '/', method, '/', bulk_name, "_fine_annot")
res_path_coarse <- paste0(res_base_path, '/', method, '/', bulk_name, "_coarse_annot")

dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)
dir.create(res_path_coarse, recursive = TRUE, showWarnings = TRUE)
dir.create(res_path_fine, recursive = TRUE, showWarnings = TRUE)

subset_list <- subset_cells(sc_matrix, sc_celltype_annotations, sc_batch, 1000, 22, 
                            sc_celltype_annotations_coarse, 
                            sc_celltype_annotations_fine)

sc_matrix <- subset_list$data
sc_celltype_annotations <- subset_list$annotations
sc_batch <- subset_list$batch_id
sc_celltype_annotations_coarse <- subset_list$annotations_coarse
sc_celltype_annotations_fine <- subset_list$annotations_fine

# Normal deconv
##############################################################################################
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

for(r in 1:replicates){

  cur_res_path <- file.path(res_path_normal, paste0('replicate_', r))
  dir.create(cur_res_path, recursive = TRUE, showWarnings = TRUE)

  bulk_matrix <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0(bulk_name, '_', bulk_norm, '.rds')))
  bulk_matrix <- as.matrix(bulk_matrix)

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

  true_fractions <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0(bulk_name, '_normal_annot_facs.rds')))

  results_list = list(
    'deconvolution' = deconvolution, 
    'true_cell_fractions' = true_fractions
  )

  saveRDS(results_list, file=paste0(cur_res_path, "/deconvolution.rds"))

}

if(method =='autogenes' | method == 'scaden'){
  unlink(signature)
}

# Coarse deconv
###############################################################################################
signature_coarse <- signature_workflow_general(
  sc_matrix, 
  sc_celltype_annotations_coarse,
  'coarse', 
  sc_dataset, 
  sc_norm, 
  sc_batch,
  method, 
  bulk_matrix,
  bulk_name, 
  bulk_norm, 
  ncores, 
  res_path_coarse,
    baseDir=baseDir
)                                               

for(r in 1:replicates){

  cur_res_path <- file.path(res_path_coarse, paste0('replicate_', r))
  dir.create(cur_res_path, recursive = TRUE, showWarnings = TRUE)

  bulk_matrix <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0(bulk_name, '_', bulk_norm, '.rds')))
  bulk_matrix <- as.matrix(bulk_matrix)

  deconvolution <- deconvolution_workflow_general(
    sc_matrix, 
    sc_celltype_annotations_coarse, 
    'coarse', 
    sc_dataset, 
    sc_norm, 
    sc_batch, 
    signature_coarse, 
    method, 
    bulk_matrix, 
    bulk_name, 
    bulk_norm, 
    ncores, 
    res_path_coarse,
    baseDir=baseDir
  )

  true_fractions <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0(bulk_name, '_coarse_annot_facs.rds')))

  results_list = list(
    'deconvolution' = deconvolution, 
    'true_cell_fractions' = true_fractions
  )

  saveRDS(results_list, file=paste0(cur_res_path, "/deconvolution.rds"))

}

if(method =='autogenes' | method == 'scaden'){
  unlink(signature_coarse)
}


# Fine deconv 
###############################################################################################

signature_fine <- signature_workflow_general(
  sc_matrix, 
  sc_celltype_annotations_fine,
  'fine', 
  sc_dataset, 
  sc_norm, 
  sc_batch,
  method, 
  bulk_matrix,
  bulk_name, 
  bulk_norm, 
  ncores, 
  res_path_fine,
  baseDir=baseDir
) 

for(r in 1:replicates){

  cur_res_path <- file.path(res_path_fine, paste0('replicate_', r))
  dir.create(cur_res_path, recursive = TRUE, showWarnings = TRUE)

  bulk_matrix <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0(bulk_name, '_', bulk_norm, '.rds')))
  bulk_matrix <- as.matrix(bulk_matrix)

  deconvolution <- deconvolution_workflow_general(
    sc_matrix, 
    sc_celltype_annotations_fine, 
    'fine', 
    sc_dataset, 
    sc_norm, 
    sc_batch, 
    signature_fine, 
    method, 
    bulk_matrix, 
    bulk_name, 
    bulk_norm, 
    ncores, 
    res_path_fine,
    baseDir=baseDir
  )

  true_fractions <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0(bulk_name, '_fine_annot_facs.rds')))

  results_list = list(
    'deconvolution' = deconvolution, 
    'true_cell_fractions' = true_fractions
  )

  saveRDS(results_list, file=paste0(cur_res_path, "/deconvolution.rds"))

}

if(method =='autogenes' | method == 'scaden'){
  unlink(signature_fine)
}
