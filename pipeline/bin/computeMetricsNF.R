#!/usr/bin/Rscript

print("Starting metric calculation ...")

library(tidyverse)


"Usage:
  computeMetricsNF.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <replicate> <subset_value> <results_dir> <baseDir>
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<bulk_name> name of simulated bulk RNAseq dataset
<bulk_path> path to simulated bulk datasets
<deconv_method>  deconv method
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<replicate> value of replicate number
<results_dir> results (base) directory
<baseDir> nextflow base directory" -> doc


args <- docopt::docopt(doc)

# store basic parameters
ncores <- as.numeric(args$ncores)
sc_dataset <- args$sc_name
sc_path <- args$sc_path
bulk_name <- args$bulk_name
bulk_path <- args$bulk_path
method <- args$deconv_method
res_base_path <- args$results_dir
subset_value <- as.numeric(args$subset_value)
replicate <- as.numeric(args$replicate)
baseDir <- args$baseDir

source(paste0(baseDir, '/bin/utils.R'))
method_normalizations <- read.table(paste0(baseDir, '/optimal_normalizations.csv'), sep = ',', header = TRUE)
sc_norm <- method_normalizations[method_normalizations$method == method, 2]
bulk_norm <- method_normalizations[method_normalizations$method == method, 3]
print(paste0('Method: ', method, '; sc-norm: ', sc_norm, '; bulk-norm: ', bulk_norm))

# load ground truth data
facs_data <- readRDS(file.path(bulk_path, bulk_name, paste0(bulk_name, '_facs.rds')))

res_path <- paste0(res_base_path, '/', method, "_", sc_dataset, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", subset_value, "_rep", replicate)

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
