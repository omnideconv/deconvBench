#!/usr/bin/Rscript

print("Starting analysis script [impact cell resolution] ...")

library(docopt)
library(Biobase)
library(omnideconv)
library(tidyverse)
reticulate::use_condaenv(condaenv = "r-omnideconv", required = TRUE)

"Usage:
  analysisNF_impact_cell_resolution.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <cell_types_fine> <replicates> <results_dir> <ncores> <baseDir> [<method_settings>]
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
<baseDir> nextflow base directory
<method_settings> settings of the method, e.g. marker_fraction=0.1 (optional)" -> doc

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
method_settings <- parse_method_settings(args$method_settings)
normalization <- method_normalization(baseDir, method)
sc_norm <- normalization$sc_norm
bulk_norm <- normalization$bulk_norm

sc_matrix <- read_sc_matrix(sc_path, sc_dataset, sc_norm)

annotation_path <- file.path(sc_path, sc_dataset, 'celltype_annotations_normal.rds')
sc_celltype_annotations <- readRDS(annotation_path)

sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))

sc_celltype_annotations_fine <- readRDS(file.path(gsub('celltype_annotations_normal.rds', 'celltype_annotations_fine.rds', annotation_path)))

coarse_file <- file.path(sc_path, sc_dataset, 'celltype_annotations_coarse.rds')


# Here we need to filter for those cell types that are in the simulated dataset. 
# NOTE: in the resolution analysis the cell types are specified in terms of finest cell types

cell_types_fine <- parse_list_argument(args$cell_types_fine, trim = TRUE)

position_vector <- sc_celltype_annotations_fine %in% cell_types_fine
sc_matrix <- sc_matrix[, position_vector]
sc_batch <- sc_batch[position_vector]
sc_celltype_annotations_fine <- sc_celltype_annotations_fine[position_vector]
sc_celltype_annotations <- sc_celltype_annotations[position_vector]


res_path_normal <- paste0(res_base_path, '/', method, '/', bulk_name, "_normal_annot")
res_path_fine <- paste0(res_base_path, '/', method, '/', bulk_name, "_fine_annot")

dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)
dir.create(res_path_fine, recursive = TRUE, showWarnings = TRUE)

if (file.exists(coarse_file)) {
  sc_celltype_annotations_coarse <- readRDS(coarse_file)
  sc_celltype_annotations_coarse <- sc_celltype_annotations_coarse[position_vector]

  res_path_coarse <- paste0(res_base_path, '/', method, '/', bulk_name, "_coarse_annot")
  dir.create(res_path_coarse, recursive = TRUE, showWarnings = TRUE)


  run_coarse <- TRUE
} else {
  run_coarse <- FALSE
}


if(run_coarse){
  subset_list <- subset_cells(sc_matrix, sc_celltype_annotations, sc_batch, 500, 22, 
                              coarse_annotations = sc_celltype_annotations_coarse, 
                              fine_annotations = sc_celltype_annotations_fine)
  
  sc_celltype_annotations_coarse <- subset_list$annotations_coarse                            

} else {
  subset_list <- subset_cells(sc_matrix, sc_celltype_annotations, sc_batch, 500, 22,  
                              fine_annotations = sc_celltype_annotations_fine)
}

sc_matrix <- subset_list$data
sc_celltype_annotations <- subset_list$annotations
sc_batch <- subset_list$batch_id
sc_celltype_annotations_fine <- subset_list$annotations_fine

# Normal deconv
##############################################################################################

# Build the signature with the reference annotated at one level and deconvolve every replicate.
# signature_bulk is only used by methods that build their signature with a bulk matrix (Scaden).
deconvolve_at_resolution <- function(annotations, level, res_path, signature_bulk){

  # The hierarchical MuSiC groups the cell types by the coarsest of the coarser levels. At the
  # coarsest level itself every cell type is its own group, which is the same as plain MuSiC.
  celltype_groups <- NULL
  if(method == 'music-hc'){
    coarser_levels <- list()
    if(level != 'coarse' && run_coarse){
      coarser_levels$coarse <- sc_celltype_annotations_coarse
    }
    if(level == 'fine'){
      coarser_levels$normal <- sc_celltype_annotations
    }
    for(coarser in coarser_levels){
      if(is.null(celltype_groups)){
        celltype_groups <- celltype_group_map(annotations, coarser)
      }
    }
    if(is.null(celltype_groups)){
      celltype_groups <- setNames(unique(annotations), unique(annotations))
    }
  }

  signature <- signature_workflow_general(
    sc_matrix, 
    annotations,
    level, 
    sc_dataset, 
    sc_norm, 
    sc_batch,
    method, 
    signature_bulk,
    bulk_name, 
    bulk_norm, 
    ncores, 
    res_path,
    baseDir=baseDir
  )

  for(r in 1:replicates){

    cur_res_path <- file.path(res_path, paste0('replicate_', r))
    dir.create(cur_res_path, recursive = TRUE, showWarnings = TRUE)

    bulk_matrix <- read_bulk_matrix(file.path(bulk_path, paste0('replicate_', r)), bulk_name, bulk_norm)

    reference <- sc_matrix
    if(method == 'rectangle'){
      # do not incluce anndata creation in runtime measurement
      reference <- rectangle_reference(sc_matrix, gsub(' ','_',annotations))
    }

    deconvolution <- deconvolution_workflow_general(
      reference, 
      annotations, 
      level, 
      sc_dataset, 
      sc_norm, 
      sc_batch, 
      signature, 
      method, 
      bulk_matrix, 
      bulk_name, 
      bulk_norm, 
      ncores, 
      res_path,
      baseDir=baseDir,
      celltype_groups = celltype_groups,
      method_settings = method_settings
    )

    true_fractions <- readRDS(file.path(bulk_path, paste0('replicate_', r), paste0(bulk_name, '_', level, '_annot_facs.rds')))

    results_list = list(
      'deconvolution' = deconvolution, 
      'true_cell_fractions' = true_fractions
    )

    saveRDS(results_list, file=paste0(cur_res_path, "/deconvolution.rds"))
  }

  if(method =='autogenes' | method == 'scaden'){
    unlink(signature)
  }
}

# The signature at normal resolution is built with the bulk matrix of the first replicate, the
# coarse and fine ones with that of the last replicate. This is how the script has always run;
# only Scaden uses the bulk matrix for its signature.
bulk_first_replicate <- read_bulk_matrix(file.path(bulk_path, 'replicate_1'), bulk_name, bulk_norm)
bulk_last_replicate <- read_bulk_matrix(file.path(bulk_path, paste0('replicate_', replicates)), bulk_name, bulk_norm)

deconvolve_at_resolution(sc_celltype_annotations, 'normal', res_path_normal, bulk_first_replicate)

if(run_coarse){
  deconvolve_at_resolution(sc_celltype_annotations_coarse, 'coarse', res_path_coarse, bulk_last_replicate)
}

deconvolve_at_resolution(sc_celltype_annotations_fine, 'fine', res_path_fine, bulk_last_replicate)
