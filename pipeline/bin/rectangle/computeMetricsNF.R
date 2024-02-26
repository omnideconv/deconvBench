#!/usr/local/bin/Rscript

library(tidyverse)
library(docopt)

print("Started metric calculation ...")
"Usage:
  computeMetricsNF.R <sc_name> <sc_norm> <bulk_dir>  <bulk_name> <bulk_norm> <deconv_method> <replicate> <subset_value> <results_dir>
Options:
<sc_name> name of sc datasets
<sc_norm> count type of sc dataset
<bulk_dir> path to bulk datasets
<bulk_name> name of bulk RNAseq dataset
<bulk_norm> normalization of RNAseq
<deconv_method>  deconv method
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<replicate> value of replicate number
<results_dir> results (base) directory" -> doc

print(doc)

args <- docopt::docopt(doc)
print(args)

sc_ds <- args$sc_name
sc_norm <- args$sc_norm
bulk_name <- args$bulk_name
bulk_norm <- args$bulk_norm
subset_value <- as.numeric(args$subset_value)
replicate <- as.numeric(args$replicate)

method <- args$deconv_method
res_base_path <- args$results_dir

facs_data <- readRDS(file.path(args$bulk_dir, args$bulk_name, paste0(args$bulk_name, '_facs.rds')))

res_path <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)

deconvolution <- readRDS(paste0(res_path, '/deconvolution.rds'))


if(args$bulk_name=='hoek' | args$bulk_name=='hoek-simulation'){
  deconvolution <- as.data.frame(deconvolution)
  deconvolution$`T cell` <- deconvolution$`T cells CD4 conv` + deconvolution$`T cells CD8` + deconvolution$`Tregs`
}


# Here we need to compute the correlation, RMSE and zRMSE
# We need to compute it by cell type, sample, all together

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
  
  results_df$estimated_frac[which(is.na(results_df$estimated_frac))] <- 0
  results_df$true_frac[which(is.na(results_df$true_frac))] <- 0
  
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

results_metric$cor_cell_type <- compute_metrics(facs_data, deconvolution, 'cor', 'cell_type')
results_metric$rmse_cell_type <- compute_metrics(facs_data, deconvolution, 'rmse', 'cell_type')

# this dataset has not enough cell-types to calculate sample-wise correlation
if(args$bulk_name != 'vanderbilt_lung'){
  results_metric$cor_sample <- compute_metrics(facs_data, deconvolution, 'cor', 'sample')
  results_metric$rmse_sample <- compute_metrics(facs_data, deconvolution, 'rmse', 'sample')
}

results_metric$facs_ground_truth <- facs_data
results_metric$deconv.results <- deconvolution

runtime_df_sig <- readRDS(paste0(res_path, '/runtime_signature.rds'))
runtime_df_dec <- readRDS(paste0(res_path, '/runtime_deconvolution.rds'))

unlink(paste0(res_path, '/runtime_signature.csv'))
unlink(paste0(res_path, '/runtime_deconvolution.csv'))

colnames(runtime_df_sig) <- c('method','sc_ds','cs_norm','bulk_ds','bulk_norm','ct_fraction','replicate','process','user.time','sys.time','elapsed')
colnames(runtime_df_dec) <- c('method','sc_ds','cs_norm','bulk_ds','bulk_norm','ct_fraction','replicate','process','user.time','sys.time','elapsed')
runtimes <- data.frame(rbind(runtime_df_sig, runtime_df_dec))

results_metric$runtimes <- runtimes

saveRDS(results_metric, paste0(res_path, '/results_metric.rds'))



## 
