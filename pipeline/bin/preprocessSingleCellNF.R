#!/usr/bin/Rscript

#path("${preProcess_dir}/${sc_ds}_${sc_norm}_perc${ct_fractions}_rep*_matrix_subsampled.rds")

library(docopt)
library(SingleCellExperiment)
library(parallel)

"Usage: 
  preprocessSingleCellNF.R <sc_matrix> <sc_annotation> <sc_batch> <sc_name> <sc_norm> <ct_fractions> <replicate> <preProcess_dir>
  
Options:
<sc_matrix> path to sc matrix
<sc_annotation> path to cell type annotation of sc matrix
<sc_batch> path batch info for sc matrix
<sc_name> name of sc datasets
<sc_norm> count type of sc dataset
<ct_fractions> fraction of cells that are subsampled from each celltype
<replicate> value of replicate number
<preProcess_dir> directory where subsampling datasets are stored" -> doc

args <- docopt::docopt(doc)
print(args)

# this ensures, that each time the script is executed, the same cells will be sampled
# still, each replicate will have different cells, as its a new loop (I tested this)
set.seed(123)


sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(args$sc_anno))
sc_batch <- readRDS(file.path(args$sc_batch))
sc_ds <- args$sc_name
sc_norm <- args$sc_norm

n_cell_types <- length(unique(sc_celltype_annotations))
n_cells_per_ct <- table(sc_celltype_annotations) 

ct_fractions <- as.numeric(args$ct_fractions)
replicate <- as.numeric(args$replicate)


cat(paste0("Preprocessing sc_matrix file with ct_fractions=", ct_fractions, " and replicate=", replicate, ")\n"))

cells_to_sample <- round(n_cells_per_ct * ct_fractions)+1 # number of cells which will be sampled from each CT (pseudocount of 1, to not have 0 cells)
print(cells_to_sample)
cell_ids <- unlist(lapply(names(cells_to_sample), function(query_ct){
  sample(x = which(sc_celltype_annotations == query_ct), size = cells_to_sample[query_ct], replace = FALSE)
}))

subset_matrix <- sc_matrix[,cell_ids]
subset_annot <- sc_celltype_annotations[cell_ids]
subset_batch <- sc_batch[cell_ids]

output_path <- paste0(args$preProcess_dir,'/',sc_ds,'_',sc_norm,'_perc',ct_fractions,'_rep',replicate)
#dir.create(output_dir, recursive = T, showWarnings = FALSE)
saveRDS(subset_matrix, paste0(output_path,'_matrix_subsampled.rds'))
saveRDS(subset_annot, paste0(output_path,'_celltype_annotations.rds'))
saveRDS(subset_batch, paste0(output_path,'_batch.rds'))

cat(paste0("Preprocessed sc_matrix file stored in '", output_path, "'\n"))


