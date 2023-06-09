#!/usr/local/bin/Rscript

library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)

print("Started deconvolution  ...")
"Usage:
  runDeconvolutionNF.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <results_dir> <run_preprocessing> <replicate> <species> <subset_value> <ncores> [<coarse>]
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
<species> type of species
<ncores> number of cores to use for method (if available)
<coarse> logical, if TRUE celltypes are mapped to higher level" -> doc

args <- docopt::docopt(doc)
print(args)

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

signature <- readRDS(paste0(res_path, "/signature.rds"))
if(method == 'scaden'){
  signature <- file.path(paste0(res_path, "/model"))
}
if(method == 'autogenes'){
  signature <- list.files(res_path, ".pickle", full.names = TRUE)
}


escapeCelltypesAutogenes <- function(celltype){
  return(gsub(" ", "xxxx", celltype))
}
reEscapeCelltypesAutogenes <- function(celltype){
  celltype <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "+", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9f", "(", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9g", ")", celltype)
  return(gsub("xxxx", " ", celltype))
}


library(Biobase) # necessary for music

###autogenes ONLY works this way!! i don't get why###
#if(method=="autogenes"){
#  library("devtools")
#  load_all("/nfs/home/extern/l.merotto/R/x86_64-pc-linux-gnu-library/4.1/omnideconv")
#}


if(method=="autogenes"){
  sc_celltype_annotations <- escapeCelltypesAutogenes(sc_celltype_annotations)
}

runtime <- system.time({
  if(method=="cpm"){
    deconvolution <- omnideconv::deconvolute(
      bulk_gene_expression = bulk_matrix, 
      signature = signature, 
      single_cell_object = sc_matrix, 
      batch_ids = sc_batch, 
      cell_type_annotations = sc_celltype_annotations, 
      method = method,
      no_cores = ncores
    )

  } else if (method=="dwls"){
    # DWLS can crash if there are not enough cells in the dataset or the cells are not differential enough between celltypes
    deconvolution <- NULL
    tryCatch(    
      deconvolution <<- omnideconv::deconvolute_dwls(
        bulk_gene_expression = bulk_matrix, 
        signature = signature, 
        dwls_submethod = 'DampenedWLS',
        verbose = TRUE
      ), error = function(e){
        print(e[['message']])
        d <- matrix(
          rep(1, ncol(bulk_matrix) * ncol(signature)),
          nrow = ncol(bulk_matrix),
          ncol = ncol(signature)
        )
        deconvolution <<- data.frame(d)
        colnames(deconvolution) <<- as.character(colnames(signature))
        rownames(deconvolution) <<- as.character(colnames(bulk_matrix))
      }
    )
    
  } else if(method=="cibersortx") {
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                          " 721a387e91c495174066462484674cb8")  
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0('/vol/omnideconv_input/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    #cx_input <- paste0('/nfs/home/students/adietrich/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    cx_output <- paste0('/vol/omnideconv_input/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    #cx_output <- paste0('/nfs/home/students/adietrich/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output)
    }

    deconvolution <- omnideconv::deconvolute_cibersortx(
      bulk_gene_expression = bulk_matrix, 
      signature = signature, 
      container = 'docker',
      verbose = TRUE,
      input_dir = cx_input,
      output_dir = cx_output
    )
    unlink(cx_input, recursive=TRUE)
    unlink(cx_output, recursive=TRUE)

  } else if (method == "cdseq"){
    source('/vol/omnideconv_input/benchmark/pipeline/bin/bin/CDseq.R')
    num_cell_type = c(5, 10, 15, 20, 30)

    deconvolution <- deconvolute_cdseq(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = sc_celltype_annotations, 
      batch_ids = sc_batch, 
      cell_type_number =  num_cell_type, 
      cpu_number = ncores,
      block_number = 5, 
      gene_subset_size = 1000
    )

    deconvolution <- t(deconvolution$cdseq_prop_merged)

  } else if (method == "bayesprism") {
    deconvolution <- omnideconv::deconvolute_bayesprism(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = sc_celltype_annotations, 
      n_cores = ncores,
      species = args$species
    )$theta
                
  } else if (method == "scaden") {
    scaden_tmp <- paste0('/vol/omnideconv_input/tmp/scaden_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    if(!dir.exists(paste0(scaden_tmp))){
      dir.create(scaden_tmp)
    }
    deconvolution <- omnideconv::deconvolute_scaden(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      temp_dir = scaden_tmp
    )
    unlink(scaden_tmp, recursive=TRUE)

  }else if (method == 'autogenes') {
    deconvolution <- omnideconv::deconvolute_autogenes(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      max_iter = 1000000
    )$proportions
    colnames(deconvolution) <- reEscapeCelltypesAutogenes(colnames(deconvolution))
    signature_dir <- paste0('/vol/omnideconv_input/tmp/autogenes_tmp_',sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate,"/")
    unlink(signature_dir, recursive = TRUE)
    file.remove(signature)
  }else {
    deconvolution <- omnideconv::deconvolute(
      bulk_gene_expression = bulk_matrix, 
      signature = signature, 
      single_cell_object = sc_matrix, 
      batch_ids = sc_batch, 
      cell_type_annotations = sc_celltype_annotations, 
      method = method
    )
  }
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
