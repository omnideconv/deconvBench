#!/usr/local/bin/Rscript
print("Started simulation script ...")

library(docopt)
library(SimBu)
library(Matrix)

"Usage: 
  simulateBulkNF.R <sc_ds> <sc_dir> <simulation_n_cells> <simulation_n_samples> <list_singlecell_datasets> <replicates> <preprocess_dir> <ncores>
Options:
<sc_ds> name of sc dataset that is used for simulations
<sc_dir> path to single cell directory
<simulation_n_cells> number of cells in each pseudo-bulk
<simulation_n_samples> number of pseudo-bulk samples
<list_singlecell_datasets> list of single cell datasets that will be used
<replicates> number of replicates to generate per dataset
<preprocess_dir> preprocessing output directory where pseudo-bulks will be stored
<ncores> number of cores for parallel simulation" -> doc
print(doc)

args <- docopt::docopt(doc)
print(args)

sc_ds <- args$sc_ds
sc_dir <- args$sc_dir
ncells <- as.numeric(args$simulation_n_cells)
nsamples <- as.numeric(args$simulation_n_samples)

datasets <- gsub('\\[|]', '', args$list_singlecell_datasets)
datasets <- strsplit(datasets , ",")[[1]]
datasets <- as.character(datasets)
print(datasets)

replicates <- as.numeric(args$replicates)
ncores <- as.numeric(args$ncores)


# We need to identify the common cell types across the datasets and generate the simulated pseudobulks with only those cell types
list.annotations <- list()
for (d in datasets){
  cur.cell.annotations <- readRDS(file.path(sc_dir, d, 'celltype_annotations.rds'))
  list.annotations[[d]] <- unique(cur.cell.annotations)
}

common.cells <- Reduce(intersect, list.annotations)
saveRDS(common.cells, file = file.path(args$preprocess_dir, 'cell_types.rds'))



for (d in datasets){
  sc_dir_cur <- paste0(sc_dir, '/', d, '/')
  sc_celltype_annotations <- readRDS(paste0(sc_dir_cur,'/celltype_annotations.rds'))

  position_vector <- sc_celltype_annotations %in% common.cells
  sc_celltype_annotations <- sc_celltype_annotations[position_vector]
  sc_matrix_raw <- readRDS(paste0(sc_dir_cur,'/matrix_counts.rds'))[, position_vector] 
  sc_matrix_norm <- readRDS(paste0(sc_dir_cur,'/matrix_norm_counts.rds'))[, position_vector] 
  
  sc_batch <- readRDS(paste0(sc_dir_cur,'/batch.rds'))[position_vector]

  simbu_ds <- SimBu::dataset(
    annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations), 
    count_matrix = sc_matrix_raw,
    tpm_matrix = sc_matrix_norm,
    name = sc_ds
  )

  for(r in 1:replicates){

    output_dir <- file.path(args$preprocess_dir, 'pseudo_bulk_impact_technology', d, paste0('/replicate_', r))

    if(dir.exists(output_dir)){
    # check if all files are present
    if(all(c('pseudo_bulk.rds','true_fractions.rds') %in% list.files(output_dir))){
        cat('Simulation with given parameters has already been done and will be skipped.')
        quit(save='no')
    }
    }else{
      dir.create(output_dir, recursive=TRUE)
    }

    simulated_bulk <- SimBu::simulate_bulk(
    data =  simbu_ds,
    scenario = 'mirror_db',
    scaling_factor = 'expressed_genes',
    total_read_counts = 1e+07,
    nsamples = nsamples,
    ncells = ncells,
    BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
    balance_even_mirror_scenario = 0.05,
    run_parallel = TRUE)
    
    saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(output_dir, '/simulation_counts.rds'))
    saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(output_dir,'/simulation_tpm.rds'))
    saveRDS(t(simulated_bulk$cell_fractions), paste0(output_dir,'/simulation_facs.rds'))
  }

}
