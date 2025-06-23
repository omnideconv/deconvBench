##### Create rds files for deconvBench #####

library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(matrixStats)


## Hao

hao_h5 <- zellkonverter::readH5AD('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao/anndata_annotated_filtered.h5ad')
celltype_annotations_original <- hao_h5$cell_type_hao
celltype_annotations_manual <- hao_h5$cell_type
matching_cells <- celltype_annotations_original == celltype_annotations_manual
hao_h5_clean <- hao_h5[,matching_cells]

# remove unepxressed genes
counts_per_gene <- Matrix::rowSums(SummarizedExperiment::assays(hao_h5_clean)[['counts']])
unexpressed_genes <- which(counts_per_gene == 0)

dataset_dir <- '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned/'
saveRDS(as.matrix(SummarizedExperiment::assays(hao_h5_clean)[['counts']][-unexpressed_genes,]), paste0(dataset_dir, '/matrix_counts.rds'))
# CPM counts
saveRDS(as.matrix(SummarizedExperiment::assays(hao_h5_clean)[['cpm']][-unexpressed_genes,]), paste0(dataset_dir, '/matrix_norm_counts.rds'))
# cell type metadata
saveRDS(hao_h5_clean$cell_type, paste0(dataset_dir, '/celltype_annotations.rds'))
# batch
saveRDS(hao_h5_clean$batch, paste0(dataset_dir, '/batch.rds'))

# SimBu dataset for simulations 
hao_ds_cleaned <- SimBu::dataset(
  annotation = data.frame(ID = colnames(hao_h5_clean), cell_type = hao_h5_clean$cell_type), 
  count_matrix = SummarizedExperiment::assays(hao_h5_clean)[['counts']][-unexpressed_genes,],
  tpm_matrix = SummarizedExperiment::assays(hao_h5_clean)[['cpm']][-unexpressed_genes,],
  name = 'haoCleaned', 
  filter_genes = F)

saveRDS(hao_ds_cleaned, paste0(dataset_dir, '/simbu_ds.rds'))

## Other datasets

basepath <- '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/'
datasets <- c("Arunachalam","BioResourceCollaboration","DominguezConde","Heimlich","Lee","SchulteSchrepping","TabulaSapiens")
res <- lapply(datasets, function(i){
  message(paste0('Handling ', i,' ...'))
  dataset_dir <- paste0(basepath, i)
  
  h5_file <- paste0(dataset_dir, '/anndata_annotated_filtered.h5ad')
  if(file.exists(h5_file)){
    ad <- zellkonverter::readH5AD(h5_file)
    message(paste0('Finished reading ', i,'.'))
    
    # raw counts
    saveRDS(as.matrix(assays(ad)[['counts']]), paste0(dataset_dir, '/matrix_counts.rds'))
    
    # CPM counts
    saveRDS(as.matrix(assays(ad)[['cpm']]), paste0(dataset_dir, '/matrix_norm_counts.rds'))
    
    # cell type metadata
    saveRDS(ad$cell_type, paste0(dataset_dir, '/celltype_annotations.rds'))
    
    # batch
    saveRDS(ad$batch, paste0(dataset_dir, '/batch.rds'))
    
    message(paste0('Finished saving ', i,', cleaning up.'))
    
    rm(ad)
    gc()  
  }
})

## Simbu ds for other datasets

basepath <- '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/'
datasets <- c("Arunachalam","BioResourceCollaboration","DominguezConde","Heimlich","Lee","SchulteSchrepping","TabulaSapiens")
res <- lapply(datasets, function(i){
  message(paste0('Handling ', i,' ...'))
  dataset_dir <- paste0(basepath, i)

  h5_file <- paste0(dataset_dir, '/anndata_annotated_filtered.h5ad')
  if(file.exists(h5_file)){
    ad <- zellkonverter::readH5AD(h5_file)
    message(paste0('Finished reading ', i,'.'))
    
    counts_per_gene <- Matrix::rowSums(SummarizedExperiment::assays(ad)[['counts']])
    unexpressed_genes <- which(counts_per_gene == 0)
    
    simbu_ds <- SimBu::dataset(
      annotation = data.frame(ID = colnames(ad), cell_type = ad$cell_type), 
      count_matrix = SummarizedExperiment::assays(ad)[['counts']][-unexpressed_genes,],
      tpm_matrix = SummarizedExperiment::assays(ad)[['cpm']][-unexpressed_genes,],
      name = i, 
      filter_genes = F)
    
    saveRDS(simbu_ds, paste0(dataset_dir, '/simbu_ds.rds'))
    
    message(paste0('Finished saving ', i,', cleaning up.'))
    
    rm(ad)
    gc()  
  }
})

