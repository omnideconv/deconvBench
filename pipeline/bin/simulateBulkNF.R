#!/usr/bin/Rscript

print("Starting simulation script ...")

library(docopt)
library(SimBu)
library(Matrix)

"Usage:
  simulateBulkNF.R <sc_ds> <sc_dir> <simulation_n_cells> <simulation_n_samples> <simulation_scenario> <simulation_bias_type> <data_dir_bulk> <ncores> [<cell_types>] [<unknown_cell_type>]
Options:
<sc_ds> name of sc dataset that is used for simulations
<sc_dir> path to single cell directory
<simulation_n_cells> number of cells in each pseudo-bulk
<simulation_n_samples> number of pseudo-bulk samples
<simulation_scenario> simulation scenario
<simulation_bias_type> type of mRNA bias to be used for simulation 
<data_dir_bulk> output dir where pseudo-bulks will be stored
<ncores> number of cores for parallel simulation)" -> doc

args <- docopt::docopt(doc)

sc_ds <- args$sc_ds
scenario <- args$simulation_scenario
bias <- args$simulation_bias_type
ncells <- as.numeric(args$simulation_n_cells)
nsamples <- as.numeric(args$simulation_n_samples)

pseudobulk_name <- paste0(sc_ds, "-ncells", ncells, "-nsamples", nsamples, "-", scenario, "-", bias, "-simulation")
output_dir <- paste0(args$data_dir_bulk, '/', pseudobulk_name)

if(dir.exists(output_dir)){
  # check if all files are present
  if(all(c(paste0(pseudobulk_name, '_counts.rds'),
           paste0(pseudobulk_name, '_tpm.rds'),
           paste0(pseudobulk_name, '_facs.rds')) 
     %in% list.files(output_dir))){

    cat('Simulation with given parameters has already been done and will be skipped.')
    quit(save='no')
  }
}else{
  dir.create(output_dir)
}

sc_dir <- paste0(args$sc_dir, '/', sc_ds, '/')
simbu_ds <- readRDS(paste0(sc_dir,'/simbu_ds.rds')) 
#sc_matrix_raw <- readRDS(paste0(sc_dir,'/matrix_counts.rds')) 
#sc_matrix_norm <- readRDS(paste0(sc_dir,'/matrix_norm_counts.rds')) 
#sc_celltype_annotations <- readRDS(paste0(sc_dir,'/celltype_annotations.rds'))
#sc_batch <- readRDS(paste0(sc_dir,'/batch.rds'))

ncores <- as.numeric(args$ncores)

#simbu_ds <- SimBu::dataset(
#  annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations), 
#  count_matrix = sc_matrix_raw,
#  tpm_matrix = sc_matrix_norm,
#  name = sc_ds, 
#  filter_genes = FALSE,
#)

set.seed(1234)
simulated_bulk <- SimBu::simulate_bulk(
  data =  simbu_ds,
  scenario = scenario,
  scaling_factor = bias,
  nsamples = nsamples,
  ncells = ncells,
  total_read_counts = 100000000,  
  BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
  run_parallel = TRUE,  
  seed = 1234
)

saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(output_dir,'/', pseudobulk_name, '_counts.rds'))
saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(output_dir,'/', pseudobulk_name, '_tpm.rds'))
saveRDS(t(simulated_bulk$cell_fractions), paste0(output_dir,'/', pseudobulk_name, '_facs.rds'))
