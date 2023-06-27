#!/usr/local/bin/Rscript

print("Started analysis for spillover script ...")

library(docopt)
library(Biobase)
library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
source('./general_functions/deconvolution_workflow_for_simulation.R')

#Sys.setenv("LD_LIBRARY_PATH"="/nfs/home/extern/l.merotto/.conda/envs/benchmark_env/x86_64-conda-linux-gnu/lib")

"Usage:
  analysisNF_spillover.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <cell_types> <results_dir> <run_preprocessing> <ncores>
Options:
<sc_matrix> path to sc matrix
<sc_annotation> path to cell type annotation of sc matrix
<sc_batch> path batch info for sc matrix
<sc_name> name of sc datasets
<sc_norm> count type of sc dataset
<bulk_dir> path to bulk datasets
<bulk_name> name of bulk RNAseq dataset
<bulk_norm> normalization of RNAseq
<deconv_method>  deconv method
<cell_types> cell types used in the simulation
<results_dir> results (base) directory
<run_preprocessing> if pre-processing has been done
<ncores> number of cores to use for method (if available)" -> doc

args <- docopt::docopt(doc)

ncores <- as.numeric(args$ncores) # in case a method can use multiple cores

sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(args$sc_anno))
sc_batch <- readRDS(file.path(args$sc_batch))
sc_ds <- args$sc_name
sc_norm <- args$sc_norm

# Here we need to filter for those cell types that are in the simulated dataset. 
# NOTE: in the resolution analysis the cell types are specified in terms of finest cell types
position_vector <- sc_celltype_annotations_fine %in% args$cell_types
sc_matrix <- sc_matrix[, position_vector]
sc_batch <- sc_batch[position_vector]
sc_celltype_annotations <- sc_celltype_annotations[position_vector]

bulk_name <- args$bulk_name
bulk_norm <- args$bulk_norm
bulk_matrix <- readRDS(file.path(args$bulk_dir, args$bulk_name, paste0(args$bulk_name, '_', args$bulk_norm, '.rds')))
bulk_matrix <- as.matrix(bulk_matrix)

method <- args$deconv_method
res_base_path <- args$results_dir

if(args$run_preprocessing == 'true'){
  subset_value <- as.numeric(args$subset_value)
  replicate <- as.numeric(args$replicate)
}else{
  subset_value <- 0
  replicate <- 0
}

res_path_normal <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm)


dir.create(res_path, recursive = TRUE, showWarnings = TRUE)


escapeCelltypesAutogenes <- function(celltype){
  celltype <- gsub("\\+", "21b2c6e87f8711ec9bf265fb9bf6ab9c", celltype)
  celltype <- gsub("-", "21b2c7567f8711ec9bf265fb9bf6ab9a", celltype)
  celltype <- gsub("\\(", "21b2c7567f8711ec9bf265fb9bf6ab9f", celltype)
  celltype <- gsub(")", "21b2c7567f8711ec9bf265fb9bf6ab9g", celltype)
  return(gsub(" ", "xxxx", celltype))
}
reEscapeCelltypesAutogenes <- function(celltype){ 
  return(gsub("xxxx", " ", celltype))
}


# Signature building 

signature <- signature_workflow_general(sc_matrix, sc_celltype_annotations, 
                                        'normal', sc_ds, sc_norm, sc_batch, bulk_matrix, method, 
                                        bulk_name, bulk_norm, res_path_normal)
# Deconvolution

deconvolution <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations, 
                                                'normal', sc_ds, sc_norm, sc_batch, signature, 
                                                method, bulk_matrix, 
                                                bulk_name, bulk_norm, res_path_normal)

saveRDS(deconvolution, file=paste0(res_path_normal, "/deconvolution_spillover_", current.cell, ".rds"))




