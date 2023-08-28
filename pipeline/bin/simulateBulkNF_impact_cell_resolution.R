#!/usr/local/bin/Rscript

print("Started simulation script ...")

library(docopt)
library(SimBu)
library(Matrix)

"Usage:
  simulateBulkNF.R <sc_ds> <sc_dir> <simulation_n_cells> <simulation_n_samples> <cell_types_minor> <preprocess_dir> <ncores> 
Options:
<sc_ds> name of sc dataset that is used for simulations
<sc_dir> path to single cell directory
<simulation_n_cells> number of cells in each pseudo-bulk
<simulation_n_samples> number of pseudo-bulk samples
<cell_types_minor> vector, subset of cell types to use for the simulation expressed in the FINER cell type annotation
<preprocess_dir> preprocessing output directory where pseudo-bulks will be stored
<ncores> number of cores for parallel simulation" -> doc

print(doc)

args <- docopt::docopt(doc)
print(args)

sc_ds <- args$sc_ds
ncells <- as.numeric(args$simulation_n_cells)
nsamples <- as.numeric(args$simulation_n_samples)
fraction_unknown <- as.numeric(args$fraction_unknown_cells)



cell_types_simulation <- gsub('\\[|]', '', args$cell_types_minor)
cell_types_simulation <- strsplit(cell_types_simulation, ",")[[1]]
print(cell_types_simulation)

pseudobulk_name <- paste0(sc_ds, '_resolution_analysis')
output_dir <- paste0(args$preprocess_dir, '/pseudo_bulk_resolution/', pseudobulk_name)

if(dir.exists(output_dir)){
  # check if all files are present
  if(all(c('pseudo_bulk.rds','true_fractions.rds') %in% list.files(output_dir))){
    cat('Simulation with given parameters has already been done and will be skipped.')
    quit(save='no')
  }
}else{
  dir.create(output_dir)
}

sc_dir <- paste0(args$sc_dir, sc_ds, '/')
sc_matrix_raw <- readRDS(paste0(sc_dir,'/matrix_counts.rds')) 
sc_matrix_norm <- readRDS(paste0(sc_dir,'/matrix_norm_counts.rds')) 
sc_batch <- readRDS(paste0(sc_dir,'/batch.rds'))

# Here we have the three level of annotations for fine, normal and coarse
sc_celltype_annotations <- readRDS(paste0(sc_dir,'celltype_annotations.rds'))
sc_celltype_annotations_fine <- readRDS(paste0(sc_dir,'celltype_annotations_fine.rds'))
sc_celltype_annotations_coarse <- readRDS(paste0(sc_dir,'celltype_annotations_coarse.rds'))


ncores <- as.numeric(args$ncores)

simbu_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(sc_matrix_raw), cell_type = sc_celltype_annotations_fine), 
  count_matrix = sc_matrix_raw,
  tpm_matrix = sc_matrix_norm,
  name = sc_ds
)

simulated_bulk <- SimBu::simulate_bulk(
  data =  simbu_ds,
  whitelist = cell_types_simulation,
  scenario = 'random',
  scaling_factor = 'expressed_genes',
  nsamples = nsamples,
  ncells = ncells,
  BPPARAM = BiocParallel::MulticoreParam(workers = ncores),
  run_parallel = TRUE
)





saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(output_dir,'/', pseudobulk_name, '_counts.rds'))
saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(output_dir,'/', pseudobulk_name, '_tpm.rds'))
saveRDS(t(simulated_bulk$cell_fractions), paste0(output_dir,'/', pseudobulk_name, 'fine_annot_facs.rds'))


# Now we need to extract the facs and group the estimates that belong to the cell types of the same families

reannotate_facs <- function(facs.table, annotation, new.annotation.level){
  
  annotation <- annotation[which(annotation$Fine %in% colnames(facs.table)), ]
  cell_types <- unique(annotation[[new.annotation.level]])
  for(c in cell_types){
    cur.cell.types <- annotation[which(annotation[[new.annotation.level]] == c), 1]
    
    if(length(cur.cell.types) > 1){
      facs.table[c] <- colSums(facs.table[, cur.cell.types])
      facs.table[, cur.cell.types] <- NULL
    }
  }
  facs.table
}

# We will use a CSV containing the information for this
table.annotations <- read.table(paste0(sc_dir, '/cell_type_mappings.csv'), header = T, sep=',')

facs.normal.annotations <- reannotate_facs(t(simulated_bulk$cell_fractions), table.annotations, 'Normal')
facs.coarse.annotations <- reannotate_facs(t(simulated_bulk$cell_fractions), table.annotations, 'Coarse')

saveRDS(facs.normal.annotations, paste0(output_dir,'/', pseudobulk_name, 'normal_annot_facs.rds'))
saveRDS(facs.coarse.annotations, paste0(output_dir,'/', pseudobulk_name, 'coarse_annot_facs.rds'))


