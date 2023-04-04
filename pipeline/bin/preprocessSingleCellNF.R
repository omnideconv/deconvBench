#!/usr/bin/Rscript

#path("${preProcess_dir}/${sc_ds}_${sc_norm}_perc${ct_fractions}_rep*_matrix_subsampled.rds")

library(docopt)
library(SingleCellExperiment)
library(parallel)

"Usage: 
  preprocessSingleCellNF.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <subset_value> <replicate> <preProcess_dir>
  
Options:
<sc_matrix> path to sc matrix
<sc_annotation> path to cell type annotation of sc matrix
<sc_batch> path batch info for sc matrix
<sc_name> name of sc datasets
<sc_norm> count type of sc dataset
<subset_value> if < 1: fraction of cell type; if > 1: number of cells per cell type
<replicate> value of replicate number
<preProcess_dir> directory where subsampling datasets are stored" -> doc

args <- docopt::docopt(doc)
print(args)

sc_ds <- args$sc_name
sc_norm <- args$sc_norm

subset_value <- as.numeric(args$subset_value)
replicate <- as.numeric(args$replicate)

output_dir <- paste0(args$preProcess_dir,'/',sc_ds,'_',sc_norm,'_perc',subset_value,'_rep',replicate)

# check if pre-processed files already exist to save time
if(dir.exists(output_dir)){
  # check if all files are present
  if(all(c('batch.rds','matrix_subsampled.rds','celltype_annotations.rds') %in% list.files(output_dir))){
    cat('Preprocessing with given parameters has already been done and will be skipped.')
    quit(save='no')
  }
}

sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(args$sc_anno))
sc_batch <- readRDS(file.path(args$sc_batch))
n_cell_types <- length(unique(sc_celltype_annotations))
n_cells_per_ct <- table(sc_celltype_annotations) 

cat(paste0("Preprocessing sc_matrix file with subset_value=", subset_value, " and replicate=", replicate, "\n"))

if(subset_value < 1){
  cells_to_sample <- round(n_cells_per_ct * subset_value)+1 # number of cells which will be sampled from each CT (pseudocount of 1, to not have 0 cells)
}else{
  cells_to_sample <- rep(subset_value, length(unique(sc_celltype_annotations)))
  names(cells_to_sample) <- names(n_cells_per_ct)
}
print(cells_to_sample)
cell_ids <- unlist(lapply(names(cells_to_sample), function(query_ct){
    if(table(sc_celltype_annotations)[query_ct] < cells_to_sample[query_ct]) {
        which(sc_celltype_annotations == query_ct) # take all cells if there are less than the subset value
    } else {
        sample(x = which(sc_celltype_annotations == query_ct), size = cells_to_sample[query_ct], replace = FALSE)
    }    
}))

subset_matrix <- sc_matrix[,cell_ids]
subset_annot <- sc_celltype_annotations[cell_ids]
subset_batch <- sc_batch[cell_ids]

dir.create(output_dir, recursive = T, showWarnings = TRUE)
print(output_dir)
saveRDS(subset_matrix, paste0(output_dir,'/matrix_subsampled.rds'))
saveRDS(subset_annot, paste0(output_dir,'/celltype_annotations.rds'))
saveRDS(subset_batch, paste0(output_dir,'/batch.rds'))

cat(paste0("Preprocessed sc_matrix file stored in '", output_dir, "'\n"))


