#!/usr/bin/Rscript

print("Starting simulation script [impact cell resolution] ...")

library(docopt)
library(SimBu)
library(Matrix)

"Usage:
  simulateBulkNF_impact_cell_resolution.R <sc_ds> <sc_dir> <simulation_n_cells> <simulation_n_samples> <cell_types_minor> <replicates> <preprocess_dir> <ncores> 
Options:
<sc_ds> name of sc dataset that is used for simulations
<sc_dir> path to single cell directory
<simulation_n_cells> number of cells in each pseudo-bulk
<simulation_n_samples> number of pseudo-bulk samples
<cell_types_minor> vector, subset of cell types to use for the simulation expressed in the FINER cell type annotation
<replicates> number of replicates to generate per dataset
<preprocess_dir> preprocessing output directory where pseudo-bulks will be stored
<ncores> number of cores for parallel simulation" -> doc


args <- docopt::docopt(doc)

sc_ds <- args$sc_ds
ncells <- as.numeric(args$simulation_n_cells)
nsamples <- as.numeric(args$simulation_n_samples)
fraction_unknown <- as.numeric(args$fraction_unknown_cells)

cell_types_simulation <- gsub('\\[|]', '', args$cell_types_minor)
cell_types_simulation <- strsplit(cell_types_simulation, ",")[[1]]
print(cell_types_simulation)

pseudobulk_name <- paste0(sc_ds, '_resolution_analysis_sim')
output_dir <- paste0(args$preprocess_dir, '/pseudo_bulk_resolution_newMethod/', pseudobulk_name)


sc_dir <- paste0(args$sc_dir, sc_ds, '/')
sc_matrix_raw <- readRDS(paste0(sc_dir,'/matrix_counts.rds')) 
sc_matrix_norm <- readRDS(paste0(sc_dir,'/matrix_norm_counts.rds')) 
sc_batch <- readRDS(paste0(sc_dir,'/batch.rds'))

# Here we have the three level of annotations for fine, normal and coarse
sc_celltype_annotations_normal <- readRDS(paste0(sc_dir,'celltype_annotations_normal.rds'))
sc_celltype_annotations_fine <- readRDS(paste0(sc_dir,'celltype_annotations_fine.rds'))
sc_celltype_annotations_coarse <- readRDS(paste0(sc_dir,'celltype_annotations_coarse.rds'))

replicates <- as.numeric(args$replicates)
ncores <- as.numeric(args$ncores)

simbu_ds_fine <- SimBu::dataset(
  annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations_fine), 
  count_matrix = sc_matrix_raw,
  tpm_matrix = sc_matrix_norm,
  name = sc_ds
)

simbu_ds_coarse <- SimBu::dataset(
  annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations_coarse), 
  count_matrix = sc_matrix_raw,
  tpm_matrix = sc_matrix_norm,
  name = sc_ds
)

simbu_ds_normal <- SimBu::dataset(
  annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations_normal), 
  count_matrix = sc_matrix_raw,
  tpm_matrix = sc_matrix_norm,
  name = sc_ds
)

# This function and the excel file will be needed for the reannotation of the ground truth fractions
reannotate_facs <- function(facs.table, annotation, new.annotation.level){
  
  facs.table <- as.data.frame(facs.table)
  annotation <- annotation[which(annotation$Fine %in% colnames(facs.table)), ]
  cell_types <- unique(annotation[[new.annotation.level]])
  for(c in cell_types){
    # These are the fine cell types
    cur.cell.types <- annotation[which(annotation[[new.annotation.level]] == c), 1]
    
    if(length(cur.cell.types) > 1){
      facs.table[c] <- rowSums(facs.table[, c(cur.cell.types)])
      facs.table[, c(cur.cell.types)] <- NULL
    }
  }
  facs.table
}

reannotate_cell_types <- function(celltypes_fine, annotation, new.annotation.level){
  
  annotation <- annotation[which(annotation$Fine %in% celltypes_fine), ]
  cell_types <- unique(annotation[[new.annotation.level]])

  cell_types
}

table.annotations <- read.table(paste0(sc_dir, '/cell_type_mappings.csv'), header = T, sep=',')


for(r in 1:replicates){

  cur_output_dir <- file.path(output_dir, paste0('/replicate_', r))

  if(dir.exists(cur_output_dir)){
    # check if all files are present
    if(all(c('pseudo_bulk.rds','true_fractions.rds') %in% list.files(output_dir))){
        cat('Simulation with given parameters has already been done and will be skipped.')
        quit(save='no')
    }
  }else{
    dir.create(cur_output_dir, recursive=TRUE)
  }

  set.seed(22+4*r)

  simulated_bulk <- SimBu::simulate_bulk(
    data =  simbu_ds_fine,
    whitelist = cell_types_simulation,
    scenario = 'mirror_db',
    scaling_factor = 'expressed_genes',
    nsamples = nsamples,
    ncells = ncells,
    balance_even_mirror_scenario = 0.05,
    BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
    run_parallel = TRUE
  )

  simulated_bulk_coarse <- SimBu::simulate_bulk(
    data =  simbu_ds_coarse,
    whitelist = reannotate_cell_types(cell_types_simulation, table.annotations, 'Coarse'),
    scenario = 'mirror_db',
    scaling_factor = 'expressed_genes',
    nsamples = nsamples,
    ncells = ncells,
    balance_even_mirror_scenario = 0.05,
    BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
    run_parallel = TRUE
  )

  simulated_bulk_normal <- SimBu::simulate_bulk(
    data =  simbu_ds_normal,
    whitelist = reannotate_cell_types(cell_types_simulation, table.annotations, 'Normal'),
    scenario = 'mirror_db',
    scaling_factor = 'expressed_genes',
    nsamples = nsamples,
    ncells = ncells,
    balance_even_mirror_scenario = 0.05,
    BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
    run_parallel = TRUE
  )

  saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(cur_output_dir,'/', pseudobulk_name, '_counts.rds'))
  saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(cur_output_dir,'/', pseudobulk_name, '_tpm.rds'))
  saveRDS(t(simulated_bulk$cell_fractions), paste0(cur_output_dir,'/', pseudobulk_name, '_fine_annot_facs.rds'))

  facs.normal.annotations <- reannotate_facs(simulated_bulk$cell_fractions, table.annotations, 'Normal')
  facs.coarse.annotations <- reannotate_facs(simulated_bulk$cell_fractions, table.annotations, 'Coarse')

  saveRDS(t(facs.normal.annotations), paste0(cur_output_dir,'/', pseudobulk_name, '_normal_annot_facs.rds'))
  saveRDS(t(facs.coarse.annotations), paste0(cur_output_dir,'/', pseudobulk_name, '_coarse_annot_facs.rds'))

}
