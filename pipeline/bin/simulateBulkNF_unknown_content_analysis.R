#!/usr/bin/Rscript

print("Starting simulation script [unkown content] ...")

library(SimBu)
library(Matrix)

"Usage:
  simulateBulkNF_unknown_content_analysis.R <sc_ds> <sc_dir> <simulation_n_cells> <simulation_n_samples> <fraction_unknown_cells> <cell_types> <unknown_cell_type> <replicate> <preprocess_dir> <ncores> 
Options:
<sc_ds> name of sc dataset that is used for simulations
<sc_dir> path to single cell directory
<simulation_n_cells> number of cells in each pseudo-bulk
<simulation_n_samples> number of pseudo-bulk samples
<fraction_unknown_cells> the fraction of unknown cellular content
<cell_types> vector, subset of cell types to use for the simulation
<unknown_cell_type> string, which cell type will be treated as unknown (i.e. not in the signature)
<replicate> number of replicate pseudobulks
<preprocess_dir> preprocessing output directory where pseudo-bulks will be stored
<ncores> number of cores for parallel simulation" -> doc


args <- docopt::docopt(doc)

sc_ds <- args$sc_ds
ncells <- as.numeric(args$simulation_n_cells)
nsamples <- as.numeric(args$simulation_n_samples)

fractions_unknown <- gsub('\\[|]', '', args$fraction_unknown_cells)
fractions_unknown <- strsplit(fractions_unknown, ",")[[1]]
fractions_unknown <- as.numeric(fractions_unknown)
print(fractions_unknown)

replicates <- as.numeric(args$replicate)

cell_types_simulation <- gsub('\\[|]', '', args$cell_types)
cell_types_simulation <- strsplit(cell_types_simulation, ",")[[1]]
print(cell_types_simulation)

unknown_cell_type <- args$unknown_cell_type
cell_types_simulation <- c(cell_types_simulation, unknown_cell_type)

pseudobulk_name <- paste0(sc_ds, '_unknown_content_sim')
output_dir <- paste0(args$preprocess_dir, '/pseudo_bulk_unknown_content/', pseudobulk_name)


if(dir.exists(output_dir)){
  # check if all files are present
  if(all(c('pseudo_bulk.rds','true_fractions.rds') %in% list.files(output_dir))){
    cat('Simulation with given parameters has already been done and will be skipped.')
    quit(save='no')
  }
}else{
  dir.create(output_dir,  recursive=TRUE)
}

sc_dir <- paste0(args$sc_dir, '/', sc_ds, '/')
sc_matrix_raw <- readRDS(paste0(sc_dir,'/matrix_counts.rds')) 
sc_matrix_norm <- readRDS(paste0(sc_dir,'/matrix_norm_counts.rds')) 
sc_celltype_annotations <- readRDS(paste0(sc_dir,'/celltype_annotations.rds'))
sc_batch <- readRDS(paste0(sc_dir,'/batch.rds'))

ncores <- as.numeric(args$ncores)

simbu_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations), 
  count_matrix = sc_matrix_raw,
  tpm_matrix = sc_matrix_norm,
  name = sc_ds
)

for (cur_cell_fraction in fractions_unknown){
  simulation_list = list()
  for (r in 1:replicates){    

    simulated_bulk <- SimBu::simulate_bulk(
      data =  simbu_ds,
      whitelist = cell_types_simulation,
      scenario = 'weighted',
      weighted_cell_type = unknown_cell_type,
      weighted_amount = cur_cell_fraction,
      scaling_factor = 'expressed_genes',
      nsamples = nsamples,
      ncells = ncells,
      BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
      run_parallel = TRUE)

    simulation_list[[as.character(cur_cell_fraction)]] <- simulated_bulk

    dir.create(paste0(output_dir, '/replicate_', r))
    saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(output_dir,'/replicate_', r, '/', pseudobulk_name, '_', cur_cell_fraction, '_counts.rds'))
    saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(output_dir,'/replicate_', r, '/', pseudobulk_name, '_', cur_cell_fraction, '_tpm.rds'))
    saveRDS(t(simulated_bulk$cell_fractions), paste0(output_dir,'/replicate_', r, '/', pseudobulk_name,  '_', cur_cell_fraction, '_facs.rds'))
  }
}
