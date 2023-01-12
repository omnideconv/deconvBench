#!/usr/bin/Rscript

print("Started signature building script ...")
library(conflicted)
conflicted::conflict_scout()
library(docopt)

#Sys.setenv("LD_LIBRARY_PATH"="/nfs/home/extern/l.merotto/.conda/envs/benchmark_env/x86_64-conda-linux-gnu/lib")

"Usage:
  computeSignatureNF.R <sc_matrix> <sc_meta> <sc_path> <rna_path> <rna_datasetname> <rna_norm> <deconv_method> <results_dir> [<coarse>]
Options:
<sc_matrix> path to sc matrix
<sc_meta> meta information of sc dataset
<sc_path> path to sc dataset
<rna_path> path to rnaseq dataset
<rna_datasetname> name of rnaseq dataset
<rna_norm> normalization of RNAseq
<deconv_method>  deconv method
<results_dir> results (base) directory
<coarse> logical, if TRUE celltypes are mapped to higher level" -> doc

n_cores <- 24 # in case a method can use multiple cores

args <- docopt::docopt(doc)

sc_path <- args$sc_path
sc_meta <- strsplit(gsub("\\]", "", gsub("\\[", "", args$sc_meta)), split = ", ")[[1]]
sc_ds <- strsplit(sc_meta[1], split = ":")[[1]][2]
sc_type <- strsplit(sc_meta[2], split = ":")[[1]][2]

sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
#sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))

rnaseq_path <- args$rna_path
rnaseq_ds <- args$rna_datasetname
rnaseq_type <- args$rna_norm
rnaseq_data <- readRDS(file.path(rnaseq_path, rnaseq_ds, paste(rnaseq_ds, "_", rnaseq_type, ".rds", sep="")))

method <- args$deconv_method
#remapping_sheet <- args$remapping_sheet
coarse <- ifelse(is.null(args$coarse), FALSE, as.logical(args$coarse))
coarse <- FALSE
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




#source("/vol/spool/bin/remapCelltypesNF.R")
#if(!coarse){
#  # regular
#  sc_dataset <- sc_ds
#  if(startsWith(sc_dataset, 'hao-')){sc_dataset = 'hao'}
#  sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
#                                                    celltype_annotations = sc_celltype_annotations, 
#                                                    method_ds = sc_dataset)
#} else {

#  # this is only designed for hoek!
#  
#}

library(Biobase)

if(method=="autogenes"){
  sc_celltype_annotations <- escapeCelltypesAutogenes(sc_celltype_annotations)
}


runtime <- system.time({
  if (method == "dwls") {
    signature <- omnideconv::build_model_dwls(
      single_cell_object = sc_matrix,
      cell_type_annotations = sc_celltype_annotations,
      dwls_method = "mast_optimized",
      ncores = n_cores
    )
  } else if (method == "cibersortx") {
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                           " 721a387e91c495174066462484674cb8") 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0('/vol/omnideconv/tmp/cibersortx_input_',sc_ds, "_", sc_type, "_", rnaseq_ds, "_", rnaseq_type)
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    cx_output <- paste0('/vol/omnideconv/tmp/cibersortx_output_',sc_ds, "_", sc_type, "_", rnaseq_ds, "_", rnaseq_type)
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
    scaden_tmp <- paste0('/vol/omnideconv/tmp/scaden_tmp_',sc_ds, "_", sc_type, "_", rnaseq_ds, "_", rnaseq_type)
    unlink(scaden_tmp, recursive=TRUE)
    if(!dir.exists(paste0(scaden_tmp))){
      dir.create(scaden_tmp)
    }
    signature <- omnideconv::build_model_scaden(
      sc_matrix,
      sc_celltype_annotations,
      rnaseq_data,
      temp_dir = scaden_tmp,
      verbose = TRUE
    )

  } else {
    signature <- omnideconv::build_model(
      single_cell_object = sc_matrix,
      cell_type_annotations = sc_celltype_annotations,
      bulk_gene_expression = rnaseq_data,
      markers = NULL,
      batch_ids = sc_batch,
      verbose = TRUE,
      method = method
    )
  }
})

res_base_path <- args$results_dir
res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_type, "_", rnaseq_ds, "_", rnaseq_type)
dir.create(res_path, recursive = TRUE, showWarnings = FALSE)


if (method %in% c("autogenes", "scaden")) { 
  #We copy the signature to the results directory
  file.copy(signature, res_path, recursive = TRUE)
  signature <- list.files(res_path, ".pickle", full.names = TRUE) 
  #rename signature
  #signature <- file.path(getwd(), basename(signature))
}

saveRDS(signature, file=paste0(res_path, "/signature.rds"), compress = FALSE)
print(runtime)
