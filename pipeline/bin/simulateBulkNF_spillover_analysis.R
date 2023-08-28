#!/usr/local/bin/Rscript

print("Started simulation script ...")

library(SimBu)
library(Matrix)
reticulate::use_miniconda(condaenv = "r-omnideconv", required = TRUE)

"Usage:
  simulateBulkNF_spillover_analysis.R <sc_ds> <sc_dir> <simulation_n_cells> <simulation_n_samples> <cell_types> <preprocess_dir> <ncores> 
Options:
<sc_ds> name of sc dataset that is used for simulations
<sc_dir> path to single cell directory
<simulation_n_cells> number of cells in each pseudo-bulk
<simulation_n_samples> number of pseudo-bulk samples
<cell_types> vector, subset of cell types to use for the simulation
<preprocess_dir> preprocessing output directory where pseudo-bulks will be stored
<ncores> number of cores for parallel simulation" -> doc

#print(doc)

args <- docopt::docopt(doc)
print(args)

sc_ds <- args$sc_ds
ncells <- as.numeric(args$simulation_n_cells)
nsamples <- as.numeric(args$simulation_n_samples)

cell_types_simulation <- gsub('\\[|]', '', args$cell_types)
cell_types_simulation <- strsplit(cell_types_simulation, ",")[[1]]
print(cell_types_simulation)

pseudobulk_name <- paste0(sc_ds, '_spillover')
output_dir <- paste0(args$preprocess_dir, '/pseudo_bulk_spillover/', pseudobulk_name, '_analysis')



if(dir.exists(output_dir)){
  # check if all files are present
  if(all(c('pseudo_bulk.rds','true_fractions.rds') %in% list.files(output_dir))){
    cat('Simulation with given parameters has already been done and will be skipped.')
    quit(save='no')
  }
}else{
  dir.create(output_dir)
}

print('Dir created')

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

print('Dataset created')

simulation_list = list()
for (cur_cell_type in cell_types_simulation){
  
  simulated_bulk <-  SimBu::simulate_bulk(
    data =  simbu_ds,
    scenario = 'pure',
    pure_cell_type = cur_cell_type,
    scaling_factor = 'expressed_genes',
    nsamples = nsamples,
    ncells = ncells,
    BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
    run_parallel = TRUE
  )

  cur_cell_type <- gsub(' ', '_', cur_cell_type)

  simulation_list[[cur_cell_type]] <- simulated_bulk
  
  saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(output_dir,'/', pseudobulk_name, '_', cur_cell_type, '_counts.rds'))
  saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(output_dir,'/', pseudobulk_name, '_', cur_cell_type, '_tpm.rds'))
  saveRDS(t(simulated_bulk$cell_fractions), paste0(output_dir,'/', pseudobulk_name, '_', cur_cell_type, '_facs.rds'))
  
}

simulated_bulk <- merge_simulations(simulation_list)

saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(output_dir,'/', pseudobulk_name, '_counts.rds'))
saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(output_dir,'/', pseudobulk_name, '_tpm.rds'))
saveRDS(t(simulated_bulk$cell_fractions), paste0(output_dir,'/', pseudobulk_name, '_facs.rds'))


