# library(tidyverse)
# library(data.table)
# 
# 
# df <- fread('revision/potential_cells.csv')
# df$collection_doi_label[which(df$collection_doi_label == 'Cambridge Institute of Therapeutic Immunology and Infectious Disease-National Institute of Health Research (CITIID-NIHR) COVID-19 BioResource Collaboration et al. (2021) Nat Med')] <- 'COVID-19 BioResource Collaboration et al. (2021) Nat Med'
# 
# df_per_doi <- df |> 
#   select(count_celltype, cell_type_broad, dataset_id_custom) |>
#   group_by(cell_type_broad, dataset_id_custom) |>
#   summarise(count_celltype = sum(count_celltype), .groups = 'drop') |>
#   group_by(dataset_id_custom) |>
#   mutate(proportion_celltype = count_celltype/sum(count_celltype)*100)
# 
# 
# ggplot(df_per_doi, aes(x=cell_type_broad, y=dataset_id_custom, fill=proportion_celltype, label=paste0(round(proportion_celltype, 3), '%')))+
#   geom_tile(color = "black")+
#   geom_text()+
#   scale_fill_gradientn(
#     colors = c("white", "lightblue", "blue")  # Adjust this color gradient as needed
#   )+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle=60, hjust = 1))
# 

##### Create rds files for deconvBench #####

library(zellkonverter)
library(SingleCellExperiment)
library(Matrix)
library(matrixStats)

#basepath <- '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/'
#dataset_dirs <- list.dirs(basepath)[-1]

# res <- lapply(dataset_dirs, function(i){
#   dataset <- basename(i)
#   dataset_dir <- i
#   
#   if(length(strsplit(dataset, '_')[[1]]) == 1){
#     message(paste0('Handling ', dataset,' ...'))
#     
#     ad <- zellkonverter::readH5AD(paste0(dataset_dir, '/anndata_annotated.h5ad'))
#     message(paste0('Finished reading ', dataset,'.'))
#     
#     # raw counts
#     saveRDS(as.matrix(assays(ad)[['counts']]), paste0(dataset_dir, '/matrix_counts.rds'))
#     
#     # CPM counts
#     saveRDS(as.matrix(assays(ad)[['cpm']]), paste0(dataset_dir, '/matrix_norm_counts.rds'))
#     
#     # cell type metadata
#     saveRDS(ad$cell_type, paste0(dataset_dir, '/celltype_annotations.rds'))
#     
#     # batch
#     saveRDS(ad$batch, paste0(dataset_dir, '/batch.rds'))
#     
#     message(paste0('Finished saving ', dataset,', cleaning up.'))
#     
#     rm(ad)
#     gc()  
#   }
# })

basepath <- '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/'
datasets <- list.files(basepath)
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


