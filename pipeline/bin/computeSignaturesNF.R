#!/usr/bin/Rscript

print("Starting signature building script ...")

library(omnideconv)
reticulate::use_condaenv(condaenv = "r-omnideconv", required = TRUE)
sessionInfo()
reticulate::py_config()

"Usage:
  computeSignaturesNF.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <results_dir> <run_preprocessing> <replicate> <subset_value> <ncores> <baseDir> 
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<bulk_name> name of simulated bulk RNAseq dataset
<bulk_path> path to simulated bulk datasets
<deconv_method>  deconv method
<results_dir> results (base) directory
<run_preprocessing> if pre-processing has been done
<replicate> value of replicate number
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<ncores> number of cores to use for method (if available)
<baseDir> nextflow base directory" -> doc


args <- docopt::docopt(doc)

# store basic parameters
ncores <- as.numeric(args$ncores)
sc_path <- args$sc_path
bulk_name <- args$bulk_name
bulk_path <- args$bulk_path
method <- args$deconv_method
res_base_path <- args$results_dir
baseDir <- args$baseDir

source(paste0(baseDir, '/bin/utils.R'))
method_normalizations <- read.table(paste0(baseDir, '/optimal_normalizations.csv'), sep = ',', header = TRUE)

# find method-specific normalizations for sc and bulk
if(method %in% method_normalizations$method){
    sc_norm <- method_normalizations[method_normalizations$method == method, 2]
    bulk_norm <- method_normalizations[method_normalizations$method == method, 3]
    print(paste0('Method: ', method, '; sc-norm: ', sc_norm, '; bulk-norm: ', bulk_norm))
} else {
    stop(paste0('Method ', method, ' is not currently supported. Exiting.')) 
}

# check if preprocessing has been performed
if(args$run_preprocessing == 'true'){
  sc_dataset <- args$sc_name
  subset_value <- gsub('perc','',strsplit(sc_dataset,'_')[[1]][2])
  replicate <- gsub('rep','',strsplit(sc_dataset,'_')[[1]][3])

  # create output directory
  res_path <- paste0(res_base_path, '/', method, "_", sc_dataset, "_", sc_norm, "_", bulk_name, "_", bulk_norm)
  dir.create(res_path, recursive = TRUE, showWarnings = TRUE)
}else{
  subset_value <- 0
  replicate <- 0
  sc_dataset <- args$sc_name

  # create output directory
  res_path <- paste0(res_base_path, '/', method, "_", sc_dataset, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
  dir.create(res_path, recursive = TRUE, showWarnings = TRUE)
}

# read scRNA-seq count matrix
if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_counts.rds'))
} else {
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_norm_counts.rds'))
}
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_dataset, 'celltype_annotations.rds'))
sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))

# read bulk expression matrix
bulk_matrix <- readRDS(file.path(args$bulk_path, bulk_name, paste0(bulk_name, '_', bulk_norm, '.rds')))
bulk_matrix <- as.matrix(bulk_matrix)




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
    res_path,
    baseDir=baseDir
  )

})

# remove old signatures for scaden
#if (method %in% c("scaden")) { 
#  # remove all pickle files in output directory with different name than current signature
#  old_signatures <- list.files(res_path, ".pickle", full.names = TRUE) 
#  if(length(old_signatures > 0)){
#    print('Found old signatures in output directory; they will be removed.')
#    sapply(old_signatures, file.remove)
#  }
#
#  #We copy the new signature to the results directory
#  file.copy(signature, res_path, recursive = TRUE)
#  signature <- list.files(res_path, ".pickle", full.names = TRUE) 
#}

# save signature object
saveRDS(signature, file=paste0(res_path, "/signature.rds"))

# measure runtime
runtime_text <- data.frame(method, 
                      args$sc_name, 
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
sessionInfo()
reticulate::py_config()
