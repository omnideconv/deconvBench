#!/usr/local/bin/Rscript

print("Started analysis for cell resolution impact script ...")

library(docopt)
library(Biobase)
library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')

#Sys.setenv("LD_LIBRARY_PATH"="/nfs/home/extern/l.merotto/.conda/envs/benchmark_env/x86_64-conda-linux-gnu/lib")

"Usage:
  analysisNF_impact_cell_resolution.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <cell_types> <results_dir> <run_preprocessing> <ncores>
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

print(args$sc_anno)
fine_type <- gsub('celltype_annotations.rds', 'celltype_annotations_fine.rds', args$sc_anno)
print(fine_type)
sc_celltype_annotations_fine <- readRDS(file.path(fine_type))
sc_celltype_annotations_coarse <- readRDS(file.path(gsub('celltype_annotations.rds', 'celltype_annotations_coarse.rds', args$sc_anno)))

# Here we need to filter for those cell types that are in the simulated dataset. 
# NOTE: in the resolution analysis the cell types are specified in terms of finest cell types
position_vector <- sc_celltype_annotations_fine %in% args$cell_types
sc_matrix <- sc_matrix[, position_vector]
sc_batch <- sc_batch[position_vector]
sc_celltype_annotations_fine <- sc_celltype_annotations_fine[position_vector]
sc_celltype_annotations_coarse <- sc_celltype_annotations_coarse[position_vector]
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

res_path_normal <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "normal_annot")
res_path_coarse <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "coarse_annot")
res_path_fine <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "fine_annot")

dir.create(res_path, recursive = TRUE, showWarnings = TRUE)
dir.create(res_path_coarse, recursive = TRUE, showWarnings = TRUE)
dir.create(res_path_fine, recursive = TRUE, showWarnings = TRUE)

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

subset_list <- subset_cells(sc_matrix, sc_celltype_annotations, sc_batch, 500, 22)

sc_matrix <- subset_list$data
sc_celltype_annotations <- subset_list$annotations
sc_batch <- subset_list$batch_id

# Signature building

signature <- signature_workflow_general(sc_matrix, sc_celltype_annotations, 
                                        'normal', sc_ds, sc_norm, sc_batch, method, bulk_matrix, 
                                        bulk_name, bulk_norm, ncores, res_path_normal)


signature.coarse <- signature_workflow_general(sc_matrix, sc_celltype_annotations_coarse, 
                                               'coarse', sc_ds, sc_norm, sc_batch, method, bulk_matrix,  
                                               bulk_name, bulk_norm, ncores, res_path_coarse)


signature.fine <- signature_workflow_general(sc_matrix, sc_celltype_annotations_fine, 
                                             'fine', sc_ds, sc_norm, sc_batch, method, bulk_matrix,  
                                              bulk_name, bulk_norm, ncores, res_path_fine)


# Deconvolution

deconvolution <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations, 
                                                'normal', sc_ds, sc_norm, sc_batch, signature, 
                                                method, bulk_matrix, bulk_name, bulk_norm, ncores,res_path_normal)


deconvolution.coarse <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations, 
                                                       'coarse', sc_ds, sc_norm, sc_batch, signature, 
                                                       method, bulk_matrix, bulk_name, bulk_norm, ncores,res_path_coarse)


deconvolution.fine <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations, 
                                                     'fine', sc_ds, sc_norm, sc_batch, signature, 
                                                     method, bulk_matrix, bulk_name, bulk_norm, ncores,res_path_fine)

saveRDS(deconvolution, file=paste0(res_path_normal, "/deconvolution.rds"))
saveRDS(deconvolution.coarse, file=paste0(res_path_coarse, "/deconvolution.rds"))
saveRDS(deconvolution.fine, file=paste0(res_path_fine, "/deconvolution.rds"))


if(method=='autogenes'){unlink(signature)}
if(method=='scaden'){unlink(signature)}











