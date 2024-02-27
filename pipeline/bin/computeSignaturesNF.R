#!/usr/bin/Rscript

print("Starting signature building script ...")

library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
source('/nfs/home/students/adietrich/omnideconv/benchmark/pipeline/bin/utils.R')

"Usage:
  computeSignatureNF.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <results_dir> <run_preprocessing> <replicate> <subset_value> <ncores> 
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
<results_dir> results (base) directory
<run_preprocessing> if pre-processing has been done
<replicate> value of replicate number
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<ncores> number of cores to use for method (if available)
<coarse> logical, if TRUE celltypes are mapped to higher level" -> doc

args <- docopt::docopt(doc)

ncores <- as.numeric(args$ncores)

sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(args$sc_anno))
sc_batch <- readRDS(file.path(args$sc_batch))
sc_ds <- args$sc_name
sc_norm <- args$sc_norm

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

# create output directory
res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
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


####################################
#### Perform signature building ####
####################################

runtime <- system.time({

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
    res_path
  )

})

if (method %in% c("scaden")) { 
  # remove all pickle files in output directory with different name than current signature
  old_signatures <- list.files(res_path, ".pickle", full.names = TRUE) 
  if(length(old_signatures > 0)){
    print('Found old signatures in output directory; they will be removed.')
    sapply(old_signatures, file.remove)
  }

  #We copy the new signature to the results directory
  file.copy(signature, res_path, recursive = TRUE)
  signature <- list.files(res_path, ".pickle", full.names = TRUE) 
}

saveRDS(signature, file=paste0(res_path, "/signature.rds"))

runtime_text <- data.frame(method, 
                      sc_ds, 
                      sc_norm, 
                      bulk_name, 
                      bulk_norm, 
                      subset_value, 
                      replicate, 
                      'SIGNATURE', 
                      runtime[['user.self']],
                      runtime[['sys.self']], 
                      runtime[['elapsed']])

saveRDS(runtime_text, file = paste0(res_path, "/runtime_signature.rds"))
