#!/usr/bin/Rscript

print("Starting metric calculation ...")

library(tidyverse)
source('/nfs/home/students/adietrich/omnideconv/benchmark/pipeline/bin/utils.R')


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

args <- docopt::docopt(doc)

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

# sum up fractions of T cell subtypes in this dataset
if(args$bulk_name=='hoek' | args$bulk_name=='hoek-simulation'){
  deconvolution <- as.data.frame(deconvolution)
  deconvolution$`T cell` <- deconvolution$`T cells CD4 conv` + deconvolution$`T cells CD8` + deconvolution$`Tregs`
}


#########################
#### Compute metrics ####
#########################

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
