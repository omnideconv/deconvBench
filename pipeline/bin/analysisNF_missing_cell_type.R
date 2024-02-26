#!/usr/local/bin/Rscript

print("Started analysis for missing cell type script ...")

library(Biobase)
library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')

"Usage:
  analysisNF_spillover.R <sc_name> <sc_path> <bulk_name> <bulk_path> <deconv_method> <missing_cell_type> <results_dir> <ncores>
Options:
<sc_name> name of sc datasets
<sc_path> path to sc dataset
<bulk_name> name of simulated bulk RNAseq dataset
<bulk_path> path to simulated bulk datasets
<deconv_method>  deconv method
<missing_cell_type> cell type to be removed from the signature
<results_dir> results (base) directory
<ncores> number of cores to use for method (if available)" -> doc

args <- docopt::docopt(doc)

sc_dataset <- args$sc_name
sc_path <- args$sc_path
bulk_name <- args$bulk_name
bulk_path <- args$bulk_path
method <- args$deconv_method
res_base_path <- args$results_dir
ncores <- as.numeric(args$ncores) # in case a method can use multiple cores

method_normalizations <- read.table('/vol/omnideconv_input/benchmark/pipeline/optimal_normalizations.csv', sep = ',', header = TRUE)
sc_norm <- method_normalizations[method_normalizations$method == method, 2]
bulk_norm <- method_normalizations[method_normalizations$method == method, 3]

if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_counts.rds'))
} else {
    sc_matrix <- readRDS(file.path(sc_path, sc_dataset, 'matrix_norm_counts.rds'))
}
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_dataset, 'celltype_annotations.rds'))
sc_batch <- readRDS(file.path(sc_path, sc_dataset, 'batch.rds'))

bulk_matrix <- readRDS(file.path(bulk_path, bulk_name, paste0(bulk_name, '_', bulk_norm, '.rds')))
bulk_matrix <- as.matrix(bulk_matrix)

# Here we need to filter for those cell types that are in the simulated dataset. 
# NOTE: in the resolution analysis the cell types are specified in terms of finest cell types

missing_cell_type <- args$missing_cell_type

if(missing_cell_type != '-'){
  position_vector <- sc_celltype_annotations != missing_cell_type
  sc_matrix <- sc_matrix[, position_vector]
  sc_batch <- sc_batch[position_vector]
  sc_celltype_annotations <- sc_celltype_annotations[position_vector]
} else {
  missing_cell_type = 'all-cells'
}




res_path_normal <- paste0(res_base_path, '/', method, "/", sc_dataset, '_', bulk_name, '_no_', missing_cell_type)
res_path_normal <- gsub(' ', '_', res_path_normal)
res_path_norml <- tolower(res_path_normal)
dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)

subset_list <- subset_cells(as.matrix(sc_matrix), sc_celltype_annotations, sc_batch, 500, 22)

sc_matrix <- subset_list$data
sc_celltype_annotations <- subset_list$annotations
sc_batch <- subset_list$batch_id

# Signature building 

signature <- signature_workflow_general(sc_matrix, sc_celltype_annotations, 
                                        'normal', sc_dataset, sc_norm, sc_batch, method, bulk_matrix,  
                                        bulk_name, bulk_norm, ncores, res_path_normal)

print('Signature built')
# Deconvolution

true_fractions <- readRDS(file.path(bulk_path, bulk_name, paste0(bulk_name, '_facs.rds')))

deconvolution <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations, 
                                                'normal', sc_dataset, sc_norm, sc_batch, signature, 
                                                method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path_normal)

colnames(deconvolution) <- gsub("xxxx", " ", colnames(deconvolution))

if(bulk_name=='hoek' | bulk_name=='hoek-simulation'){
  deconvolution <- as.data.frame(deconvolution)
  deconvolution$`T cell` <- switch(missing_cell_type,
    'T cells CD4 conv' = (deconvolution$`T cells CD8` + deconvolution$`Tregs`),
    'T cells CD8' = (deconvolution$`T cells CD4 conv` + deconvolution$`Tregs`),
    'Tregs' = (deconvolution$`T cells CD4 conv` + deconvolution$`T cells CD8`),
    (deconvolution$`T cells CD4 conv` + deconvolution$`T cells CD8` + deconvolution$`Tregs`)
  )
}

results_metric <- list()
results_metric$cor_cell_type <- compute_metrics(true_fractions, deconvolution, 'cor', 'cell_type')
results_metric$rmse_cell_type <- compute_metrics(true_fractions, deconvolution, 'rmse', 'cell_type')
results_metric$facs_ground_truth <- true_fractions
results_metric$deconv.results <- deconvolution

saveRDS(results_metric, file=paste0(res_path_normal, "/deconvolution.rds"))

if(method =='autogenes' | method == 'scaden'){
  unlink(signature)
}



