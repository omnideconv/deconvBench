#!/usr/local/bin/Rscript

print("Started analysis for cell resolution impact script ...")

library(docopt)
library(Biobase)
library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')

"Usage:
  analysisNF_impact_cell_resolution.R <deconv_method> <results_dir> <ncores>
Options:
<deconv_method>  deconv method
<results_dir> results (base) directory
<ncores> number of cores to use for method (if available)" -> doc

args <- docopt::docopt(doc)

method <- args$deconv_method
res_base_path <- args$results_dir
ncores <- as.numeric(args$ncores) # in case a method can use multiple cores

bulk_dataset <- 'hoek'
sc_dataset <- 'hao-sampled-1'

method_normalizations <- read.table('/vol/omnideconv_input/benchmark/pipeline/optimal_normalizations.csv', sep = ',', header = TRUE)
sc_norm <- method_normalizations[method_normalizations$method == method, 2]
bulk_norm <- method_normalizations[method_normalizations$method == method, 3]

sc_path <- file.path('/vol/omnideconv_input/omnideconv_data/singleCell/', sc_dataset)

if(sc_norm == 'counts'){
    sc_matrix <- readRDS(file.path(sc_path, 'matrix_counts.rds'))
} else {
    sc_matrix <- readRDS(file.path(sc_path, 'matrix_norm_counts.rds'))
}

sc_celltype_annotations <- readRDS(file.path(sc_path, 'celltype_annotations.rds'))
sc_celltype_annotations_coarse <- readRDS(file.path(sc_path, 'celltype_annotations_coarse.rds'))

sc_batch <- readRDS(file.path(sc_path, 'batch.rds'))

bulk_path <- file.path('/vol/omnideconv_input/omnideconv_data/PBMC', bulk_dataset)
bulk_matrix <- readRDS(file.path(bulk_path, paste0(bulk_dataset, '_', bulk_norm, '.rds')))
bulk_matrix <- as.matrix(bulk_matrix)


res_path_normal <- paste0(res_base_path, "/normal_annot/", method, '/')
res_path_coarse <- paste0(res_base_path, "/coarse_annot/", method, '/')

dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)
dir.create(res_path_coarse, recursive = TRUE, showWarnings = TRUE)

# Signature building

signature <- signature_workflow_general(sc_matrix, sc_celltype_annotations, 
                                        'normal', sc_dataset, sc_norm, sc_batch, method, bulk_matrix, 
                                        bulk_dataset, bulk_norm, ncores, res_path_normal)


deconvolution <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations, 
                                                'normal', sc_dataset, sc_norm, sc_batch, signature, 
                                                method, bulk_matrix, bulk_dataset, bulk_norm, ncores, res_path_normal)

deconvolution <- as.data.frame(deconvolution)

deconvolution$`T cells` <- deconvolution$`T cells CD4 conv` + deconvolution$`T cells CD8` + deconvolution$`Tregs` 


true_fractions <- readRDS(file.path(bulk_path, paste0(bulk_dataset, '_facs.rds')))

results_list = list(
    'deconvolution' = deconvolution, 
    'true_cell_fractions' = true_fractions
)

saveRDS(results_list, file=paste0(res_path_normal, "/deconvolution.rds"))

if(method =='autogenes' | method == 'scaden'){
  unlink(signature)
}

signature_coarse <- signature_workflow_general(sc_matrix, sc_celltype_annotations_coarse, 
                                               'coarse', sc_dataset, sc_norm, sc_batch, method, bulk_matrix,  
                                               bulk_dataset, bulk_norm, ncores, res_path_coarse)


deconvolution_coarse <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations_coarse, 
                                                       'coarse', sc_dataset, sc_norm, sc_batch, signature_coarse, 
                                                       method, bulk_matrix, bulk_dataset, bulk_norm, ncores, res_path_coarse)

results_list_coarse = list(
    'deconvolution' = deconvolution_coarse, 
    'true_cell_fractions' = true_fractions
)

saveRDS(results_list_coarse, file=paste0(res_path_coarse, "/deconvolution.rds"))

if(method =='autogenes' | method == 'scaden'){
  unlink(signature_coarse)
}
