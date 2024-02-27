#!/usr/bin/Rscript

print("Starting deconvolution script ...")

library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
source('/nfs/home/students/adietrich/omnideconv/benchmark/pipeline/bin/utils.R')

"Usage:
  runDeconvolutionNF.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <results_dir> <run_preprocessing> <replicate> <subset_value> <species> <ncores> 
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
<species> type of species
<ncores> number of cores to use for method (if available)
<coarse> logical, if TRUE celltypes are mapped to higher level" -> doc

args <- docopt::docopt(doc)

ncores <- as.numeric(args$ncores) # in case a method can use multiple cores

sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(args$sc_anno))
sc_batch <- readRDS(file.path(args$sc_batch))
sc_ds <- args$sc_name
method <- args$deconv_method

method_normalizations <- read.table('/nfs/home/students/adietrich/omnideconv/benchmark/pipeline/optimal_normalizations.csv', sep = ',', header = TRUE)
sc_norm <- method_normalizations[method_normalizations$method == method, 2]
bulk_norm <- method_normalizations[method_normalizations$method == method, 3]
#TODO: check if specified normalization is available

if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_counts.rds'))
} else {
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_norm_counts.rds'))
}

bulk_name <- args$bulk_name
bulk_matrix <- readRDS(file.path(args$bulk_dir, args$bulk_name, paste0(args$bulk_name, '_', args$bulk_norm, '.rds')))
bulk_matrix <- as.matrix(bulk_matrix)

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

# path to temporary directory, if needed
tmp_dir_path <- paste0(res_base_path, '/tmp_', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)

# load signature
signature <- readRDS(paste0(res_path, "/signature.rds"))
if(method == 'scaden'){
  signature <- file.path(paste0(res_path, "/model"))
}

# escape celltype names for python based method autogenes
if(method=="autogenes"){
  sc_celltype_annotations <- escapeCelltypesAutogenes(sc_celltype_annotations)
}

###############################
#### Perform devonvolution ####
###############################

runtime <- system.time({

  deconvolution <- deconvolution_workflow_general(
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

colnames(deconvolution) <- gsub('\\.',' ',colnames(deconvolution))
saveRDS(deconvolution, file=paste0(res_path, "/deconvolution.rds"))

runtime_text <- data.frame(method, 
                      sc_ds, 
                      sc_norm, 
                      bulk_name, 
                      bulk_norm, 
                      subset_value, 
                      replicate, 
                      'DECONVOLUTION', 
                      runtime[['user.self']], 
                      runtime[['sys.self']], 
                      runtime[['elapsed']])

saveRDS(runtime_text, file = paste0(res_path, "/runtime_deconvolution.rds"))
