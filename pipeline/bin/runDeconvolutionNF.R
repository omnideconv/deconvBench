#!/usr/bin/Rscript

print("Starting deconvolution script ...")

library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)

"Usage:
  runDeconvolutionNF.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <results_dir> <run_preprocessing> <replicate> <subset_value> <species> <ncores> <baseDir>
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
<species> type of species
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
sc_norm <- method_normalizations[method_normalizations$method == method, 2]
bulk_norm <- method_normalizations[method_normalizations$method == method, 3]
print(paste0('Method: ', method, '; sc-norm: ', sc_norm, '; bulk-norm: ', bulk_norm))

# check if preprocessing has been performed
if(args$run_preprocessing == 'true'){
  subset_value <- as.numeric(args$subset_value)
  replicate <- as.numeric(args$replicate)
  sc_dataset <- paste0(args$sc_name,'_',sc_norm,'_perc',subset_value,'_rep',replicate)
}else{
  subset_value <- 0
  replicate <- 0
  sc_dataset <- args$sc_name
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

# create output directory
res_path <- paste0(res_base_path, '/', method, "_", args$sc_name, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
dir.create(res_path, recursive = TRUE, showWarnings = TRUE)

# load signature
signature <- readRDS(paste0(res_path, "/signature.rds"))
if(method == 'scaden'){
  signature <- file.path(paste0(res_path, "/model"))
}

# escape celltype names for python based method autogenes
if(method=="autogenes"){
  sc_celltype_annotations <- escapeCelltypesAutogenes(sc_celltype_annotations)
}

# set cibersortx batch correction options
if(method=='cibersortx'){
  datasets_technologies <- read.table(paste0(baseDir, '/sc_datasets_technologies.csv'), sep = ',', header = TRUE)
  cur_tech <- datasets_technologies[datasets_technologies$sc_dataset == args$sc_name, 2]

  if(grepl("simulation" , bulk_name)){
      s_mode <- FALSE
      b_mode <- FALSE
  } else if(cur_tech!='10X'){
      s_mode <- FALSE
      b_mode <- TRUE
  } else {
      s_mode <- TRUE
      b_mode <- FALSE
  }
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
    signature,
    method, 
    bulk_matrix,
    bulk_name, 
    bulk_norm, 
    ncores, 
    res_path,
    rmbatch_S_mode = s_mode,
    rmbatch_B_mode = b_mode,
    baseDir=baseDir
  )

})

# save deconvolution result
colnames(deconvolution) <- gsub('\\.',' ',colnames(deconvolution))
saveRDS(deconvolution, file=paste0(res_path, "/deconvolution.rds"))

# measure runtime
runtime_text <- data.frame(method, 
                      args$sc_name, 
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
