#!/usr/local/bin/Rscript

print("Started simulation script ...")

library(docopt)
library(SimBu)
library(Matrix)

"Usage:
  simulateBulkNF.R <sc_ds> <sc_norm> <sc_dir> <simulation_n_cells> <simulation_n_samples> <simulation_scenario> <preprocess_dir> <ncores> 
Options:
<sc_ds> name of sc dataset that is used for simulations
<sc_norm> normalization of sc dataset that is used for simulations
<sc_dir> path to single cell directory
<simulation_n_cells> number of cells in each pseudo-bulk
<simulation_n_samples> number of pseudo-bulk samples
<simulation_scenario> simulation scenario
<preprocess_dir> preprocessing output directory where pseudo-bulks will be stored
<ncores> number of cores for parallel simulation" -> doc

print(doc)

args <- docopt::docopt(doc)
print(args)

sc_ds <- args$sc_ds
sc_norm <- args$sc_norm

scenario <- args$simulation_scenario
ncells <- args$simulation_n_cells
nsamples <- args$simulation_n_samples

output_dir <- paste0(args$preprocess_dir, '/pseudo_bulk/', sc_ds, "_", sc_norm, "_ncells", ncells, "_nsamples", nsamples, "_", scenario)

if(dir.exists(output_dir)){
  # check if all files are present
  if(all(c('pseudo_bulk.rds','true_fractions.rds') %in% list.files(output_dir))){
    cat('Simulation with given parameters has already been done and will be skipped.')
    quit(save='no')
  }
}

sc_dir <- paste0(args$sc_dir, '/', sc_ds, '/')
sc_matrix_raw <- readRDS(paste0(sc_dir,'/matrix_counts.rds')) 
sc_matrix_norm <- readRDS(paste0(sc_dir,'/matrix_norm_counts.rds')) 
sc_celltype_annotations <- readRDS(paste0(sc_dir,'/celltype_annotations.rds'))
sc_batch <- readRDS(paste0(sc_dir,'/batch.rds'))

ncores <- args$ncores

simbu_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations), 
  count_matrix = sc_matrix_raw,
  tpm_matrix = sc_matrix_norm,
  name = sc_ds
)

simulated_bulk <- SimBu::simulate_bulk(
  data =  simbu_ds,
  scenario = scenario,
  scaling_factor = 'expressed_genes',
  nsamples = nsamples,
  ncells = ncells,
  BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
  run_parallel = TRUE
)


saveRDS(simulated_bulk$bulk, output_dir)
saveRDS(simulated_bulk$cell_fractions, output_dir)


