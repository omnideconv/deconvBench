## create subsets of Hao, TabulaSapiens and BioResource datasets

library(tidyverse)


### Hao
set.seed(234)
anno_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao/celltype_annotations.rds') 
anno_hao_tbl <- anno_hao |>
  as.tibble() |>
  mutate(index = 1:length(anno_hao)) |>
  group_by(value)

anno_hao_subset <- slice_sample(anno_hao_tbl, prop = .1, replace = FALSE)

counts_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao/matrix_counts.rds') 
cpm_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao/matrix_norm_counts.rds') 
batch_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao/batch.rds') 

counts_hao_subset <- counts_hao[,anno_hao_subset$index]
cpm_hao_subset <- cpm_hao[,anno_hao_subset$index]
batch_hao_subset <- batch_hao[anno_hao_subset$index]

if(!dir.exists('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao-sampled')){
  dir.create('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao-sampled')
}

saveRDS(batch_hao_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao-sampled/batch.rds')
saveRDS(counts_hao_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao-sampled/matrix_counts.rds')
saveRDS(cpm_hao_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao-sampled/matrix_norm_counts.rds')
saveRDS(as.character(anno_hao_subset$value), '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/Hao-sampled/celltype_annotations.rds')

### TabulaSapiens
set.seed(234)
anno_TabulaSapiens <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens/celltype_annotations.rds') 
anno_TabulaSapiens_tbl <- anno_TabulaSapiens |>
  as.tibble() |>
  mutate(index = 1:length(anno_TabulaSapiens)) |>
  group_by(value)

anno_TabulaSapiens_subset <- slice_sample(anno_TabulaSapiens_tbl, prop = .1, replace = FALSE)

counts_TabulaSapiens <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens/matrix_counts.rds') 
cpm_TabulaSapiens <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens/matrix_norm_counts.rds') 
batch_TabulaSapiens <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens/batch.rds') 

counts_TabulaSapiens_subset <- counts_TabulaSapiens[,anno_TabulaSapiens_tbl$index]
cpm_TabulaSapiens_subset <- cpm_TabulaSapiens[,anno_TabulaSapiens_tbl$index]
batch_TabulaSapiens_subset <- anno_TabulaSapiens[anno_TabulaSapiens_tbl$index]

if(!dir.exists('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens-sampled')){
  dir.create('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens-sampled')
}

saveRDS(batch_TabulaSapiens_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens-sampled/batch.rds')
saveRDS(counts_TabulaSapiens_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens-sampled/matrix_counts.rds')
saveRDS(cpm_TabulaSapiens_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens-sampled/matrix_norm_counts.rds')
saveRDS(as.character(anno_TabulaSapiens_subset$value), '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/TabulaSapiens-sampled/celltype_annotations.rds')

### BioResource
set.seed(234)
anno_BioResource <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration//celltype_annotations.rds') 
anno_BioResource_tbl <- anno_BioResource |>
  as.tibble() |>
  mutate(index = 1:length(anno_BioResource)) |>
  group_by(value)

anno_BioResource_subset <- slice_sample(anno_BioResource_tbl, prop = .1, replace = FALSE)

counts_BioResource <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration/matrix_counts.rds') 
cpm_BioResource <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration/matrix_norm_counts.rds') 
batch_BioResource <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration/batch.rds') 

counts_BioResource_subset <- counts_BioResource[,anno_BioResource_tbl$index]
cpm_BioResource_subset <- cpm_BioResource[,anno_BioResource_tbl$index]
batch_BioResource_subset <- anno_BioResource[anno_BioResource_tbl$index]

if(!dir.exists('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration-sampled')){
  dir.create('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration-sampled')
}

saveRDS(batch_BioResource_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration-sampled/batch.rds')
saveRDS(counts_BioResource_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration-sampled/matrix_counts.rds')
saveRDS(cpm_BioResource_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration-sampled/matrix_norm_counts.rds')
saveRDS(as.character(anno_BioResource_subset$value), '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/BioResourceCollaboration-sampled/celltype_annotations.rds')