## create subsets of Hao and other datasets

library(tidyverse)

#### Hao ####
set.seed(234)
anno_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned/celltype_annotations.rds') 
anno_hao_tbl <- anno_hao |>
  as.tibble() |>
  mutate(index = 1:length(anno_hao)) |>
  group_by(value)

anno_hao_subset <- anno_hao_tbl |>
  group_by(value) |>
  group_modify(~ {
    n_total <- nrow(.x)
    n_sample <- max(10, ceiling(n_total * 0.1))   # 10% of each celltype
    slice_sample(.x, n = min(n_sample, n_total))  # Avoid oversampling small groups
  }) |>
  ungroup()

counts_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned/matrix_counts.rds') 
cpm_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned/matrix_norm_counts.rds') 
batch_hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned/batch.rds') 

counts_hao_subset <- counts_hao[,anno_hao_subset$index]
cpm_hao_subset <- cpm_hao[,anno_hao_subset$index]
batch_hao_subset <- batch_hao[anno_hao_subset$index]

if(!dir.exists('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled')){
  dir.create('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled')
}

saveRDS(batch_hao_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/batch.rds')
saveRDS(counts_hao_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/matrix_counts.rds')
saveRDS(cpm_hao_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/matrix_norm_counts.rds')
saveRDS(as.character(anno_hao_subset$value), '/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned-sampled/celltype_annotations.rds')

#### other datasets ####
set.seed(234)

names <- c('Arunachalam','BioResourceCollaboration','DominguezConde','Heimlich','Lee','SchulteSchrepping','TabulaSapiens')
target_total <- 15000
min_n_cells <- 20

out <- lapply(names, function(i){
  print(paste0('Handling ',i,' ...'))
  dir <- paste0('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/',i)
  dir_sampled_out <- paste0(dir,'-sampled')
  
  anno <- readRDS(paste0(dir,'/celltype_annotations.rds'))
  
  anno_tbl <- anno |>
    as.tibble() |>
    mutate(index = 1:length(anno)) |>
    group_by(value)
  
  celltype_counts <- anno_tbl %>%
    dplyr::count(value, name = "n")
  
  total_cells <- sum(celltype_counts$n)
  
  # Sample proportionally
  anno_subset <- anno_tbl %>%
    inner_join(celltype_counts, by = "value") %>%
    group_by(value) %>%
    group_modify(~ {
      n_total <- nrow(.x)
      target_n1 = round(n_total / total_cells * target_total)
      target_n2 = max(target_n1, min_n_cells)
      target_n_final = min(target_n2, n_total)
      slice_sample(.x, n = target_n_final, replace = FALSE)
    }) %>%
    ungroup()
  
  counts_subset <- readRDS(paste0(dir,'/matrix_counts.rds'))[,anno_subset$index]
  cpm_subset <- readRDS(paste0(dir,'/matrix_norm_counts.rds'))[,anno_subset$index]
  batch_subset <- readRDS(paste0(dir,'/batch.rds'))[anno_subset$index]
  
  # remove genes with 0 expression in subset
  counts_per_gene <- rowSums(counts_subset)
  unexpressed_genes <- which(counts_per_gene == 0)
  
  counts_subset_expressed_only <- counts_subset[-unexpressed_genes,]
  cpm_subset_expressed_only <- cpm_subset[-unexpressed_genes,]
  
  if(!dir.exists(dir_sampled_out)){
    dir.create(dir_sampled_out)
  }
  
  saveRDS(batch_subset, paste0(dir_sampled_out, '/batch.rds'))
  saveRDS(counts_subset_expressed_only, paste0(dir_sampled_out, '/matrix_counts.rds'))
  saveRDS(cpm_subset_expressed_only, paste0(dir_sampled_out, '/matrix_norm_counts.rds'))
  saveRDS(as.character(anno[anno_subset$index]), paste0(dir_sampled_out, '/celltype_annotations.rds'))
})

#### Myers ####

set.seed(234)
anno_myers <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/celltype_annotations.rds') 
anno_myers_tbl <- anno_myers |>
  as.tibble() |>
  mutate(index = 1:length(anno_myers)) |>
  group_by(value)

anno_myers_subset <- anno_myers_tbl |>
  group_by(value) |>
  group_modify(~ {
    n_total <- nrow(.x)
    n_sample <- max(10, ceiling(n_total * 0.1))   # 10% of each celltype
    slice_sample(.x, n = min(n_sample, n_total))  # Avoid oversampling small groups
  }) |>
  ungroup()

counts_myers <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/matrix_counts.rds') 
cpm_myers <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/matrix_norm_counts.rds') 
batch_myers <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/batch.rds') 

myers_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(counts_myers), cell_type = anno_myers), 
  count_matrix = counts_myers,
  tpm_matrix = cpm_myers,
  name = 'Myers', 
  filter_genes = T)
saveRDS(myers_ds, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/simbu_ds.rds')

counts_myers_subset <- counts_myers[,anno_myers_subset$index]
cpm_myers_subset <- cpm_myers[,anno_myers_subset$index]
batch_myers_subset <- batch_myers[anno_myers_subset$index]

if(!dir.exists('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers-sampled')){
  dir.create('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers-sampled')
}

saveRDS(batch_myers_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers-sampled/batch.rds')
saveRDS(counts_myers_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers-sampled/matrix_counts.rds')
saveRDS(cpm_myers_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers-sampled/matrix_norm_counts.rds')
saveRDS(as.character(anno_myers_subset$value), '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers-sampled/celltype_annotations.rds')


#### Allen ####

set.seed(234)
anno_allen <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen/celltype_annotations.rds') 
anno_allen_tbl <- anno_allen |>
  as.tibble() |>
  mutate(index = 1:length(anno_allen)) |>
  group_by(value)

anno_allen_subset <- anno_allen_tbl |>
  group_by(value) |>
  group_modify(~ {
    n_total <- nrow(.x)
    n_sample <- max(10, ceiling(n_total * 0.1))   # 10% of each celltype
    slice_sample(.x, n = min(n_sample, n_total))  # Avoid oversampling small groups
  }) |>
  ungroup()

counts_allen <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen/matrix_counts.rds') 
cpm_allen <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen/matrix_norm_counts.rds') 
batch_allen <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen/batch.rds') 

allen_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(counts_allen), cell_type = anno_allen), 
  count_matrix = counts_allen,
  tpm_matrix = cpm_allen,
  name = 'Allen', 
  filter_genes = T)
saveRDS(allen_ds, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen/simbu_ds.rds')

counts_allen_subset <- counts_allen[,anno_allen_subset$index]
cpm_allen_subset <- cpm_allen[,anno_allen_subset$index]
batch_allen_subset <- batch_allen[anno_allen_subset$index]

if(!dir.exists('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen-sampled')){
  dir.create('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen-sampled')
}

saveRDS(batch_allen_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen-sampled/batch.rds')
saveRDS(counts_allen_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen-sampled/matrix_counts.rds')
saveRDS(cpm_allen_subset, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen-sampled/matrix_norm_counts.rds')
saveRDS(as.character(anno_allen_subset$value), '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/allen-sampled/celltype_annotations.rds')
