library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
library(cowplot)

#source('./cell_palette.R')
#source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')
methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')

process_results_df <- function(res, method, resolution, replicate, dataset, predicted = TRUE){
  
  res <- res %>%
    tibble::rownames_to_column(., var='sample') %>%
    gather(., key='celltype', value = 'predicted_value', -'sample') %>%
    #mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., method = method) %>%
    mutate(., resolution = resolution) %>%
    mutate(., replicate = replicate) %>% 
    mutate(., dataset = dataset)
  
  if(!predicted){
    colnames(res)[colnames(res) == 'predicted_value'] <- 'true_value' 
  }
  
  res
}

cell_palette <- c('B cells'='#999933',
                  'Macrophages'='#CC6677',
                  'mDC'='#882255',
                  'Monocytes'='#AA4499',
                  'NK cells'='#DDCC77',
                  'T cell'='#332288',
                  'T cells CD4 conv'='#117733',
                  'T cells CD8'='#44AA99',
                  'Tregs'='#88CCEE',
                  'pDC'='#8D4B00',
                  'Neutrophils'='#FFFB6A',
                  'Stromal cells' = '#288BA8',
                  'Macrophages-Monocytes'='#CC6677',
                  'T and NK cells'='#332288',
                  'DCs'='#8D4B00')

# 1: List directories, methods, cell types
#############################################################################
resolution.deconv.results <- list.files('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_resolution', full.names=F, recursive=T)

resolution.deconv.results.lamb <- resolution.deconv.results[grep('lambrechts', resolution.deconv.results)]
resolution.deconv.results.lamb <- resolution.deconv.results.lamb[grep('deconvolution.rds', resolution.deconv.results.lamb)]

resolution.deconv.results.wu <- resolution.deconv.results[grep('wu', resolution.deconv.results)]
resolution.deconv.results.wu <- resolution.deconv.results.wu[grep('deconvolution.rds', resolution.deconv.results.wu)]

resolution.deconv.results <- c(resolution.deconv.results.lamb, resolution.deconv.results.wu)

metadata.table <- resolution.deconv.results %>%
  tibble(path = ., 
         method = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         dataset_level = map_vec(., function(x) strsplit(x, split = '/')[[1]][2]),
         replicate = map_vec(., function(x) strsplit(x, split = '/')[[1]][3])) %>%
  mutate(dataset_level = gsub("_resolution_analysis_sim|_annot", "", dataset_level),
         replicate = gsub('replicate_', '', replicate)) %>%
  separate(dataset_level, c("dataset", "level"), sep="_")

#2: Combine these in a unique dataframe
################################################################################

#table.annotations <- read.table(paste0('/vol/omnideconv_input/omnideconv_data/singleCell/lambrechts', '/cell_type_mappings.csv'), header = T, sep=',')
#colnames(table.annotations) <- c('fine','normal','coarse')

tmp <- lapply(1:nrow(metadata.table), function(i){
  result <- readRDS(paste0('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_resolution/', metadata.table$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  colnames(result) <- gsub("xxxx", " ", colnames(result))
  colnames(result) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(result))
  colnames(result) <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "-", colnames(result))
  
  true.fractions <- readRDS(paste0('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_resolution/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("xxxx", " ", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "-", colnames(true.fractions))
  
  result <- result %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i], metadata.table$replicate[i], metadata.table$dataset[i])
  
  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i], 
                       metadata.table$replicate[i], metadata.table$dataset[i], predicted = FALSE)
  
  left_join(result, true.fractions)
  
})

data <- data.table::rbindlist(tmp)

################################################################################
# Heatmaps for correlation, RMSE
################################################################################

correlation.results <- data %>%
  group_by(celltype, method, resolution, dataset) %>%
  dplyr::summarize(corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate, 
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results.all <- data %>%
  group_by(method, resolution, dataset) %>%
  dplyr::summarize(celltype = 'all',
                   corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate, 
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results <- rbind(correlation.results, correlation.results.all)

correlation.results$resolution <- factor(correlation.results$resolution, levels=c('coarse', 'normal', 'fine'))
correlation.results$celltype <- factor(correlation.results$celltype, 
                                       levels = c(unique(correlation.results$celltype)))

ggplot(correlation.results, aes(x=celltype, y=method, fill=corr))+
  geom_tile()+
  geom_text(aes(label = round(corr,2)), size = 2.5)+
  facet_grid(.~resolution, scales = "free_x", space = "free_x")+
  theme_minimal()+
  scale_fill_gradient2(low = '#ce273f', high='#2d87bb', limits=c(-1,1))+
  rotate_x_text(angle=60)+
  coord_flip()

ggplot(correlation.results, aes(x=celltype, y=method, fill=rmse)) +
  geom_tile()+
  geom_text(aes(label = round(rmse,2)), size = 2.5)+
  facet_grid(.~resolution, scales = "free_x", space = "free_x")+
  theme_minimal()+
  scale_fill_gradient(low='white',high = 'purple') +
  rotate_x_text(angle=60)

