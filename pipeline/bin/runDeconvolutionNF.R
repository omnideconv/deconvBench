#!/usr/bin/Rscript

print("Starting deconvolution script ...")

library(omnideconv)
reticulate::use_condaenv(condaenv = "r-omnideconv", required = TRUE)
sessionInfo()
reticulate::py_config()

"Usage:
  runDeconvolutionNF.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <results_dir> <run_preprocessing> <replicate> <subset_value> <species> <ncores> <baseDir> <signature_dir> [<method_settings>]
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<bulk_name> name of simulated bulk RNAseq dataset
<bulk_path> path to simulated bulk datasets
<deconv_method>  deconv method
<results_dir> directory of this result; signature and deconvolution files are read from and written to it
<run_preprocessing> if pre-processing has been done
<replicate> value of replicate number
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<species> type of species
<ncores> number of cores to use for method (if available)
<baseDir> nextflow base directory
<signature_dir> directory with the files of the signature step; does not exist if no signature step was run
<method_settings> settings of the method, e.g. marker_fraction=0.1 (optional)" -> doc

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
method_settings <- parse_method_settings(args$method_settings)
normalization <- method_normalization(baseDir, method)
sc_norm <- normalization$sc_norm
bulk_norm <- normalization$bulk_norm
print(paste0('Method: ', method, '; sc-norm: ', sc_norm, '; bulk-norm: ', bulk_norm))

# check if preprocessing has been performed
if(args$run_preprocessing == 'true'){
  sc_dataset <- args$sc_name
  subset_value <- gsub('perc','',strsplit(sc_dataset,'_')[[1]][2])
  replicate <- gsub('rep','',strsplit(sc_dataset,'_')[[1]][3])

  # create output directory
  res_path <- res_base_path
  dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
  res_path <- normalizePath(res_path)
}else{
  subset_value <- 0
  replicate <- 0
  sc_dataset <- args$sc_name

  # create output directory
  res_path <- res_base_path
  dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
  res_path <- normalizePath(res_path)
}

# read scRNA-seq count matrix
sc_matrix <- read_sc_matrix(sc_path, sc_dataset, sc_norm)
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_dataset, 'celltype_annotations.rds'))
# the names before any method-specific escaping
reference_celltypes <- unique(as.character(sc_celltype_annotations))
# groups of cell types, only used by the hierarchical MuSiC
celltype_groups <- read_celltype_groups(method, file.path(sc_path, sc_dataset), sc_celltype_annotations)
sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))

# read bulk expression matrix
bulk_matrix <- read_bulk_matrix(file.path(args$bulk_path, bulk_name), bulk_name, bulk_norm)

# take over the files of the signature step, so that the result directory is complete
if(dir.exists(args$signature_dir)){
  signature_files <- setdiff(list.files(args$signature_dir), 'tmp')
  file.copy(file.path(args$signature_dir, signature_files), res_path, recursive = TRUE)
}

if(file.exists(paste0(res_path, "/runtime_signature.rds"))){
  # a signature that is shared between bulk datasets gets the name of the bulk dataset of this result
  runtime_signature <- readRDS(paste0(res_path, "/runtime_signature.rds"))
  runtime_signature[[4]] <- bulk_name
  saveRDS(runtime_signature, file = paste0(res_path, "/runtime_signature.rds"))
} else {
  # no signature step was run, because this method does not build a signature;
  # the (empty) signature and its runtime are recorded here
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
  saveRDS(signature, file=paste0(res_path, "/signature.rds"))
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
}

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
  } else if(cur_tech != '10X'){
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

if(method == 'rectangle'){
  # do not incluce anndata creation in runtime measurement
  sc_matrix <- rectangle_reference(sc_matrix, sc_celltype_annotations)
}

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
    rmbatch_S_mode = FALSE,
    rmbatch_B_mode = FALSE,
    baseDir=baseDir,
    celltype_groups = celltype_groups,
    method_settings = method_settings
  )

})

# save deconvolution result, with the cell-type names of the reference
colnames(deconvolution) <- restore_celltype_names(colnames(deconvolution), reference_celltypes)
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
