#!/usr/local/bin/Rscript

print("Started analysis for spillover script ...")

library(docopt)
library(Biobase)
library(omnideconv)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)
source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')

"Usage:
  analysisNF_spillover.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <bulk_dir> <bulk_name> <bulk_norm> <deconv_method> <fraction_unknown_cells> <cell_types> <unknown_cell_type> <results_dir> <run_preprocessing> <ncores>
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
<fraction_unknown_cells> the fraction of unknown cellular content
<cell_types> vector, subset of cell types to use for the simulation
<unknown_cell_type> string, which cell type will be treated as unknown (i.e. not in the signature)
<results_dir> results (base) directory
<run_preprocessing> if pre-processing has been done
<ncores> number of cores to use for method (if available)" -> doc

args <- docopt::docopt(doc)
print(args)

ncores <- as.numeric(args$ncores) # in case a method can use multiple cores

sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(args$sc_anno))
sc_batch <- readRDS(file.path(args$sc_batch))
sc_ds <- args$sc_name
sc_norm <- args$sc_norm

unknown_cell_type <- args$unknown_cell_type

fractions_unknown <- gsub('\\[|]', '', args$fraction_unknown_cells)
fractions_unknown <- strsplit(fractions_unknown, ",")[[1]]
fractions_unknown <- as.numeric(fractions_unknown)
print(fractions_unknown)


# Here we need to filter for those cell types that are in the simulated dataset. 
# NOTE: in the resolution analysis the cell types are specified in terms of finest cell types

cell_types_simulation <- gsub('\\[|]', '', args$cell_types)
cell_types_simulation <- strsplit(cell_types_simulation, ",")[[1]]
print(cell_types_simulation)

position_vector <- sc_celltype_annotations %in% cell_types_simulation
sc_matrix <- sc_matrix[, position_vector]
sc_batch <- sc_batch[position_vector]
sc_celltype_annotations <- sc_celltype_annotations[position_vector]

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

res_path_normal <- paste0(res_base_path, '/', method, "_", sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm)

# Docker needs that
res_path_normal <- tolower(res_path_normal)

dir.create(res_path_normal, recursive = TRUE, showWarnings = TRUE)


subset_list <- subset_cells(sc_matrix, sc_celltype_annotations, sc_batch, 500, 22)

sc_matrix <- subset_list$data
sc_celltype_annotations <- subset_list$annotations
sc_batch <- subset_list$batch_id

# Signature building 

signature <- signature_workflow_general(sc_matrix, sc_celltype_annotations, 
                                        'normal', sc_ds, sc_norm, sc_batch, method, bulk_matrix, 
                                        bulk_name, bulk_norm, ncores,res_path_normal)
# Deconvolution

for(cur_cell_fraction in fractions_unknown){
  
  bulk_matrix <- readRDS(file.path(args$bulk_dir, args$bulk_name, paste0(args$bulk_name, '_', cur_cell_fraction, '_', args$bulk_norm, '.rds')))
  bulk_matrix <- as.matrix(bulk_matrix)
  
  deconvolution <- deconvolution_workflow_general(sc_matrix, sc_celltype_annotations, 
                                                  'normal', sc_ds, sc_norm, sc_batch, signature, 
                                                  method, bulk_matrix, 
                                                  bulk_name, bulk_norm, ncores, res_path_normal)
  
  saveRDS(deconvolution, file=paste0(res_path_normal, "/deconvolution__", unknown_cell_type, '_', cur_cell_fraction, ".rds"))
  
}

saveRDS(deconvolution, file=paste0(res_path_normal, "/deconvolution.rds"))

if(method=='autogenes'){unlink(signature)}
if(method=='scaden'){unlink(signature)}

