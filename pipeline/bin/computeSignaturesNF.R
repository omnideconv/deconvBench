#!/usr/bin/Rscript

print("Started signature building script ...")

library(docopt)
library(Biobase)
library(omnideconv)

#Sys.setenv("LD_LIBRARY_PATH"="/nfs/home/extern/l.merotto/.conda/envs/benchmark_env/x86_64-conda-linux-gnu/lib")

"Usage:
  computeSignatureNF.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <results_dir> <run_preprocessing> <replicate> <subset_value> <ncores> [<coarse>]
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
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<replicate> value of replicate number
<ncores> number of cores to use for method (if available)
<coarse> logical, if TRUE celltypes are mapped to higher level" -> doc

args <- docopt::docopt(doc)

ncores <- as.numeric(args$ncores) # in case a method can use multiple cores

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

res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
dir.create(res_path, recursive = TRUE, showWarnings = FALSE)

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

runtime <- system.time({
  if (method == "dwls") {
    signature <- omnideconv::build_model_dwls(
      single_cell_object = sc_matrix,
      cell_type_annotations = sc_celltype_annotations,
      dwls_method = "mast_optimized",
      ncores = ncores
    )
  } else if (method == "cibersortx") {
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                           " 721a387e91c495174066462484674cb8") 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0('/vol/omnideconv/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    #cx_input <- paste0('/nfs/home/students/adietrich/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    cx_output <- paste0('/vol/omnideconv/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    #cx_output <- paste0('/nfs/home/students/adietrich/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output)
    }
    signature <- omnideconv::build_model_cibersortx(
      sc_matrix, 
      sc_celltype_annotations, 
      container = "docker", 
      verbose = TRUE, 
      input_dir = cx_input, 
      output_dir = cx_output
    )
    
  } else if(method == "scaden"){
    scaden_tmp <- paste0('/vol/omnideconv/tmp/scaden_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    unlink(scaden_tmp, recursive=TRUE)
    if(!dir.exists(paste0(scaden_tmp))){
      dir.create(scaden_tmp)
    }
    signature <- omnideconv::build_model_scaden(
      sc_matrix,
      sc_celltype_annotations,
      bulk_matrix,
      temp_dir = scaden_tmp,
      verbose = TRUE
    )
  } else if(method == "autogenes"){
    sc_celltype_annotations <- escapeCelltypesAutogenes(sc_celltype_annotations)
    
    signature_dir <- paste0('/vol/omnideconv/tmp/autogenes_tmp_',sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate,"/")
    unlink(signature_dir, recursive=TRUE)
    if(!dir.exists(paste0(signature_dir))){
      dir.create(signature_dir)
    }
    signature <- omnideconv::build_model_autogenes(
      sc_matrix,
      sc_celltype_annotations,
      output_dir = signature_dir,
      verbose = TRUE
    )
  } else {
    signature <- omnideconv::build_model(
      single_cell_object = sc_matrix,
      cell_type_annotations = sc_celltype_annotations,
      bulk_gene_expression = bulk_matrix,
      markers = NULL,
      batch_ids = sc_batch,
      verbose = TRUE,
      method = method
    )
  }
})


if (method %in% c("autogenes", "scaden")) { 
  #We copy the signature to the results directory
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


