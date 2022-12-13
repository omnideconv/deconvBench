#!/usr/bin/Rscript

"Usage: 
  runDeconvolutionNF.R <sc_matrix> <sc_meta>  <sc_path> <rna_path> <rna_datasetname> <rna_norm> <deconv_method> <results_dir> [<coarse>]
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
rnaseq_data <- as.matrix(rnaseq_data)

method <- args$deconv_method

res_base_path <- args$results_dir
res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_type, "_", rnaseq_ds, "_", rnaseq_type)


signature <- readRDS(paste0(res_path, "/signature.rds"))

##remap celltype_annotations##
#source("/vol/spool/bin/remapCelltypesNF.R")

#sc_dataset <- sc_ds
#if(startsWith(sc_dataset, 'hao-sampled')){sc_dataset = 'hao'}

#sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
#                                                  celltype_annotations = sc_celltype_annotations, 
#                                                  method_ds = sc_dataset)

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
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = rnaseq_data, 
                                           signature = signature, 
                                           single_cell_object = sc_matrix, 
                                           batch_ids = sc_batch, 
                                           cell_type_annotations = sc_celltype_annotations, 
                                           method = method,
                                           no_cores = n_cores)

} else if(method=="cibersortx") {
  omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                         " 721a387e91c495174066462484674cb8")  
  # CibersortX does not work with tmp directories in a Docker in Docker setup
  # --> created fixed input and output directories!
  cx_input <- '/vol/spool/tmp/cibersortx_input'
  if(!dir.exists(paste0(cx_input))){
    dir.create(cx_input)
  }
  cx_output <- '/vol/spool/tmp/cibersortx_output'
      if(!dir.exists(paste0(cx_output))){
    dir.create(cx_output)
  }

  deconvolution <- deconvolute_cibersortx(bulk_gene_expression = rnaseq_data, 
                                          signature = signature, 
                                          container = 'docker',
                                          verbose = TRUE,
                                          input_dir = cx_input,
                                          output_dir = cx_output)

  unlink(cx_input)
  unlink(cx_output)

} else if (method == "cdseq"){
  source('/vol/omnideconv/benchmark/pipeline/bin/bin/CDseq.R')
  num_cell_type = c(5, 10, 15, 20, 30)

  deconvolution <- deconvolute_cdseq(bulk_gene_expression = rnaseq_data, 
                      single_cell_object = sc_matrix, 
                      cell_type_annotations = sc_celltype_annotations, 
                      batch_ids = sc_batch, 
                      cell_type_number =  num_cell_type, 
                      cpu_number = n_cores,
                      block_number = 5, 
                      gene_subset_size = 1000)

  deconvolution <- t(deconvolution$cdseq_prop_merged)
  #deconvolution <- omnideconv::normalize_deconv_results(deconvolution)
    

} else if (method == "bayesprism") {
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = rnaseq_data, 
              signature = signature, 
              single_cell_object = sc_matrix, 
              batch_ids = sc_batch, 
              cell_type_annotations = sc_celltype_annotations, 
              method = method,
              n_cores = n_cores)
              
} else {
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = rnaseq_data, 
              signature = signature, 
              single_cell_object = sc_matrix, 
              batch_ids = sc_batch, 
              cell_type_annotations = sc_celltype_annotations, 
              method = method)
}






})

res_base_path <- args$results_dir
res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_type, "_", rnaseq_ds, "_", rnaseq_type)



saveRDS(deconvolution, file=paste0(res_path, "/deconvolution.rds"))


#saveRDS(deconvolution, 
#        paste("deconvolution_", method, "_", sc_ds, "_", sc_type, "_", rnaseq_ds, "_", rnaseq_type, ".rds", sep=""), 
#        compress = FALSE)



# Here we need to compute the correlation, RMSE and zRMSE
# We need to compute it by cell type, sample, all together

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
  deconvolution$`T cell` <- deconvolution$`T cell CD4+` + deconvolution$`T cell CD8+` + deconvolution$`T cell regulatory (Tregs)`

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




saveRDS(results_metric,
        paste0(res_path, '/results_metric.rds'))


print(runtime)
