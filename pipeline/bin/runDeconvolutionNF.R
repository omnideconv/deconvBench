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
  if(method=="autogenes"){

    deconvolution <- omnideconv::deconvolute_autogenes( 
      single_cell_object = sc_matrix,
      cell_type_annotations = sc_celltype_annotations,
      bulk_gene_expression = bulk_matrix, 
      max_iter = 1000000,
      ngen = 5000,
      verbose = TRUE
    )$proportions
    colnames(deconvolution) <- reEscapeCelltypesAutogenes(colnames(deconvolution))

  } else if (method=="bayesprism"){
    
    # need to use different parameter in case of simulated bulk data
    if(grepl("simulation", bulk_name)){
      update_gibbs <- FALSE
      which_theta <- 'first'
    }else{
      update_gibbs <- TRUE
      which_theta <- 'final'
    }
    
    deconvolution <- omnideconv::deconvolute_bayesprism(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = sc_celltype_annotations, 
      n_cores = ncores,
      species = args$species,
      update_gibbs = update_gibbs,
      which_theta = which_theta
    )$theta
    
  } else if (method=="bisque"){
    deconvolution <- t(omnideconv::deconvolute_bisque(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = sc_celltype_annotations,
      batch_ids = sc_batch,
      verbose = TRUE,
    )$bulk.props)
    
  } else if (method=="cibersortx"){

    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                          " 721a387e91c495174066462484674cb8")  
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(tmp_dir_path,'_input')
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    cx_input <- paste0(tmp_dir_path,'_output')
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output)
    }

    # batch correction options for cibersortx
    rmbatch_B_mode <- FALSE
    rmbatch_S_mode <- TRUE
    if(grepl("simulation" , bulk_name)){
      rmbatch_S_mode <- FALSE
    }
    
    deconvolution <- omnideconv::deconvolute_cibersortx(
      single_cell_object = sc_matrix,
      bulk_gene_expression = bulk_matrix,
      cell_type_annotations = sc_celltype_annotations, 
      signature = signature,
      container = 'docker',
      verbose = TRUE,
      input_dir = cx_input,
      output_dir = cx_output,
      rmbatch_B_mode = rmbatch_B_mode,
      rmbatch_S_mode = rmbatch_S_mode
    )
    unlink(cx_input, recursive=TRUE)
    unlink(cx_output, recursive=TRUE)
    
  } else if(method=="dwls") {
    
    deconvolution <- NULL
    tryCatch(    
      deconvolution <<- omnideconv::deconvolute_dwls(
        bulk_gene_expression = bulk_matrix, 
        signature = signature, 
        dwls_submethod = 'DampenedWLS',
        verbose = TRUE
      ), error = function(e){
        # DWLS can crash if there are not enough cells in the dataset or the cells are not differential enough between celltypes
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

  } else if (method == "music"){
    deconvolution <- omnideconv::deconvolute_music(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = sc_celltype_annotations,
      batch_ids = sc_batch,
      verbose = TRUE,
    )$Est.prop.weighted

  } else if (method == "scdc") {    

                
  } else if (method == "scaden") {
    scaden_tmp <- paste0('/vol/omnideconv_results/tmp/scaden_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)
    if(!dir.exists(paste0(scaden_tmp))){
      dir.create(scaden_tmp)
    }
    deconvolution <- omnideconv::deconvolute_scaden(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      temp_dir = scaden_tmp
    )
    unlink(scaden_tmp, recursive=TRUE)

  } else {
    message('Selected method is not supported in the benchmark. Please check again.')
    stop()
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
