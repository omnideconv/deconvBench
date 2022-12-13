#!/usr/bin/Rscript

library(conflicted)
conflicted::conflict_scout()
library(docopt)


"Usage:
  computeMetricsNF.R <sc_meta> <rnaseq> <rnaseq_norm> <mode> <results_dir> <facs> 
Options:
<sc_meta> meta information of sc dataset
<rnaseq> name of rnaseq dataset
<rnaseq_norm> normalization of RNAseq
<mode> deconv method
<results_dir> results (base) directory
<facs> folder of matrix with facs fractions" -> doc

print(doc)



args <- docopt::docopt(doc)
#results <- args$deconvolution

#resultVec <- strsplit(gsub("\\[", "", gsub("\\]", "", results)), ", ")[[1]]
#print(resultVec)


# We set up the results directory
sc_dataset <- args$sc_meta
bulk_dataset <- args$rnaseq
bulk_norm <- args$rnaseq_norm
cur_method <- args$mode
facs_folder <- args$facs

facs_folder <- file.path(strsplit(gsub("\\[", "", gsub("\\]", "", facs_folder)), ", ")[[1]])


sc_meta <- strsplit(gsub("\\]", "", gsub("\\[", "", args$sc_meta)), split = ", ")[[1]]
sc_ds <- strsplit(sc_meta[1], split = ":")[[1]][2]
sc_norm <- strsplit(sc_meta[2], split = ":")[[1]][2]
#print(sc_ds)

##remap celltype annotations of facs##
#source("/vol/spool/bin/remapCelltypesNF.R")



facs_path_new <- file.path(facs_folder, bulk_dataset, paste(bulk_dataset, "_facs.rds", sep=""))
print(facs_path_new)
facs_data <- readRDS(facs_path_new)

#if(endsWith(bulk_dataset, '-simulation')){
#  bulk_dataset_correspondance <- 'hao'
#}

#celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet,
#                                               celltype_annotations = rownames(facs_data),
#                                               method_ds = bulk_dataset_correspondance)
#rownames(facs_data) <- celltype_annotations

#print(results)


res_base_path <- args$results_dir
res_path <- paste0(res_base_path, '/', cur_method, "_", sc_ds, "_", sc_norm, "_", bulk_dataset, "_", bulk_norm)

deconvolution.results <- readRDS(paste0(res_path, '/deconvolution.rds'))

print('done')

if(bulk_dataset=='hoek'){
  deconvolution.results <- as.data.frame(deconvolution.results)
  deconvolution.results$`T cell` <- deconvolution.results$`T cell CD4+` + deconvolution.results$`T cell CD8+` + deconvolution.results$`T cell regulatory (Tregs)`

}

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



results_metric <- list()

if(cur_method !='cdseq'){
  results_metric$cor_sample <- compute_metrics(facs_data, deconvolution.results, 'cor', 'sample')
  results_metric$cor_cell_type <- compute_metrics(facs_data, deconvolution.results, 'cor', 'cell_type')
  results_metric$rmse_sample <- compute_metrics(facs_data, deconvolution.results, 'rmse', 'sample')
  results_metric$rmse_cell_type <- compute_metrics(facs_data, deconvolution.results, 'rmse', 'cell_type')
}

results_metric$facs_groud_truth <- facs_data
results_metric$deconv.results <- deconvolution.results



res_base_path <- args$results_dir
res_path <- paste0(res_base_path, '/', cur_method, "_", sc_ds, "_", sc_norm, "_", bulk_dataset, "_", bulk_norm)



saveRDS(results_metric,
        paste0(res_path, '/results_metric.rds'))
