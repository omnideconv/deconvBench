#!/usr/bin/Rscript

print("Started deconvolution  ...")
"Usage:
  runDeconvolutionNF.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <results_dir> <run_preprocessing> <replicate> <ct_fractions> <ncores> [<coarse>]
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
<ct_fractions> fraction of cells that are subsampled from each celltype
<replicate> value of replicate number
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

method <- args$deconv_method
res_base_path <- args$results_dir

if(args$run_preprocessing == 'true'){
  ct_fractions <- as.numeric(args$ct_fractions)
  replicate <- as.numeric(args$replicate)
}else{
  ct_fractions <- 0
  replicate <- 0
}

res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", ct_fractions, "_rep", replicate)

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

  } else if(method=="cibersortx") {
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                          " 721a387e91c495174066462484674cb8")  
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    #cx_input <- paste0('/vol/omnideconv/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", ct_fractions, "_rep", replicate)
    cx_input <- paste0('/nfs/home/students/adietrich/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", ct_fractions, "_rep", replicate)
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    #cx_output <- paste0('/vol/omnideconv/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", ct_fractions, "_rep", replicate)
    cx_output <- paste0('/nfs/home/students/adietrich/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", ct_fractions, "_rep", replicate)
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
    source('/vol/omnideconv/benchmark/pipeline/bin/bin/CDseq.R')
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
    deconvolution <- omnideconv::deconvolute(
      bulk_gene_expression = bulk_matrix, 
      signature = signature, 
      single_cell_object = sc_matrix, 
      batch_ids = sc_batch, 
      cell_type_annotations = sc_celltype_annotations, 
      method = method,
      n_cores = ncores
    )
                
  } else if (method == "scaden") {
    scaden_tmp <- paste0('/vol/omnideconv/tmp/scaden_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", ct_fractions, "_rep", replicate)
    if(!dir.exists(paste0(scaden_tmp))){
      dir.create(scaden_tmp)
    }
    deconvolution <- omnideconv::deconvolute_scaden(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      temp_dir = scaden_tmp
    )
    unlink(scaden_tmp)
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

res_base_path <- args$results_dir
res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", ct_fractions, "_rep", replicate)

saveRDS(deconvolution, file=paste0(res_path, "/deconvolution.rds"))


# Here we need to compute the correlation, RMSE and zRMSE
# We need to compute it by cell type, sample, all together

print("Started metric calculation ...")

library(tidyverse)

compute_metrics <- function(ref_facs, deconvolution, metric = c('cor', 'rmse'),
                            comparison = c('sample', 'cell_type')){
  
  compute_rmse <- function(x, y, zscored = FALSE, weighted = FALSE){
    if(zscored == TRUE){
      sd_x = sd(x)
      sd_y = sd(y)
      if(sd_x == 0){sd_x = 1}
      if(sd_y == 0){sd_y = 1}
      x <- (x - mean(x))/(sd_x)
      y <- (y - mean(y))/(sd_y)
    }
    
    rmse <- sqrt(mean((x - y)^2)) 
    if(weighted == TRUE){rmse <- rmse/max(y)}
    return(rmse)
  }
  
  
  
  deconvolution <- as.data.frame(deconvolution) %>%
    rownames_to_column(., 'sample') %>%
    gather(., key = 'cell_type', value = 'estimated_frac', -sample)
  
  ref_facs <- as.data.frame(ref_facs) %>%
    rownames_to_column(., 'cell_type') %>%
    gather(., key = 'sample', value = 'true_frac', -cell_type)
  
  results_df <- deconvolution %>%
    inner_join(., ref_facs) #, by = c('cell_type' ,'sample'))
  
  
  if(metric == 'cor'){
    res_metric <- results_df %>%
      group_by(!!! syms(comparison)) %>%
      dplyr::summarize(cor = cor.test(estimated_frac, true_frac)$estimate,
                       pval = cor.test(estimated_frac, true_frac)$p.value)
    
    cor_all_elements <- data.frame(comp = 'all' ,
                                   'cor' = cor.test(results_df$estimated_frac, results_df$true_frac)$estimate,
                                   'pval' = cor.test(results_df$estimated_frac, results_df$true_frac)$p.value,
                                   row.names = (nrow(res_metric) +1))
    colnames(cor_all_elements)[1] <- comparison
    res_metric <- rbind(res_metric, cor_all_elements)
    
  } else if(metric == 'rmse'){
    res_metric <- results_df %>%
      group_by(!!! syms(comparison)) %>%
      dplyr::summarize(RMSE = compute_rmse(estimated_frac, true_frac),
                       zRMSE = compute_rmse(estimated_frac, true_frac, zscored=TRUE), 
                       wRMSE = compute_rmse(estimated_frac, true_frac, weighted=TRUE)
                       )
    
    rmse_all_elements <- data.frame(comp = 'all' ,
                                    'RMSE' = compute_rmse(results_df$estimated_frac, results_df$true_frac),
                                    'zRMSE' = compute_rmse(results_df$estimated_frac, results_df$true_frac, zscored=TRUE),
                                    'wRMSE' = compute_rmse(results_df$estimated_frac, results_df$true_frac, weighted=TRUE),
                                    row.names = (nrow(res_metric) +1))
    
    colnames(rmse_all_elements)[1] <- comparison
    res_metric <- rbind(res_metric, rmse_all_elements)
  }
  
  return(res_metric)
  
}


facs_folder <- file.path(strsplit(gsub("\\[", "", gsub("\\]", "", rnaseq_path)), ", ")[[1]])
facs_path_new <- file.path(facs_folder, rnaseq_ds, paste(rnaseq_ds, "_facs.rds", sep=""))
print(facs_path_new)
facs_data <- readRDS(facs_path_new)

if(rnaseq_ds=='hoek'){
  deconvolution <- as.data.frame(deconvolution)
  deconvolution$`T cell` <- deconvolution$`T cells CD4 conv` + deconvolution$`T cells CD8` + deconvolution$`Tregs`
}

results_metric <- list()

if(method !='cdseq'){
  results_metric$cor_sample <- compute_metrics(facs_data, deconvolution, 'cor', 'sample')
  results_metric$cor_cell_type <- compute_metrics(facs_data, deconvolution, 'cor', 'cell_type')
  results_metric$rmse_sample <- compute_metrics(facs_data, deconvolution, 'rmse', 'sample')
  results_metric$rmse_cell_type <- compute_metrics(facs_data, deconvolution, 'rmse', 'cell_type')
}

results_metric$facs_groud_truth <- facs_data
results_metric$deconv.results <- deconvolution

saveRDS(results_metric, paste0(res_path, '/results_metric.rds'))
print(runtime)
