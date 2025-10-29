library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
library(cowplot)

# 1: List directories, methods, cell types
#############################################################################
impact.technology.results <- list.files('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_impact_technology', full.names=F, recursive=T)

#impact.technology.results <- impact.technology.results[grep('lambrechts', impact.technology.results)]
impact.technology.results <- impact.technology.results[grep('deconvolution.rds', impact.technology.results)]

impact.technology.results <- impact.technology.results[grep('vanderbilt', impact.technology.results)]

metadata.table <- impact.technology.results %>%
  tibble(path = .,
         method_pseudobulk = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         pseudobulk_signature = map_vec(., function(x) strsplit(x, split = '/')[[1]][2])) %>%
  separate(method_pseudobulk, c("method", "dataset_pseudobulk"), sep="_") %>%
  separate(pseudobulk_signature, c("dataset", "tissue", "bulk", "sc_dataset", "signature"))

metadata.table[, c(4, 5, 6, 8)] <- NULL


#2: Combine these in a unique dataframe
################################################################################
data <- NULL

process_results_df <- function(res, method, pseudobulk_ds, signature_ds, predicted = TRUE){

  res <- res %>%
    tibble::rownames_to_column(., var='sample') %>%
    gather(., key='celltype', value = 'predicted_value', -'sample') %>%
    mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., method = method) %>%
    mutate(., dataset_pseudobulk = pseudobulk_ds) %>%
    mutate(., dataset_signature = signature_ds)

  if(!predicted){
    colnames(res)[colnames(res) == 'predicted_value'] <- 'true_value'
  }

  res
}

tmp <- lapply(1:nrow(metadata.table), function(i){
  result <- readRDS(paste0('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_impact_technology/', metadata.table$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  colnames(result) <- gsub("xxxx", " ", colnames(result))
  colnames(result) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(result))
  
  result$`T cells CD4` <- result$`T cells CD4 conv` + result$`Tregs`
  result$`T cells CD4 conv` <- NULL
  result$`Tregs` <- NULL
  result$`T cells` <- result$`T cells CD4` + result$`T cells CD8`
  
  true.fractions <- readRDS(paste0('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_impact_technology/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("xxxx", " ", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))
  
  true.fractions$`T cells CD4` <- true.fractions$`T cells CD4 conv`
  true.fractions$`T cells CD4 conv` <- NULL
  if("CD4" %in% colnames(true.fractions)){
    true.fractions$`T cells CD4` <- true.fractions$`CD4`
    true.fractions$`CD4` <- NULL
  }
  if("CD8" %in% colnames(true.fractions)){
    true.fractions$`T cells CD8` <- true.fractions$`CD8`
    true.fractions$`CD8` <- NULL
  }

  
  true.fractions$`T cells` <- true.fractions$`T cells CD4` + true.fractions$`T cells CD8`
  
  result <- result %>%
    process_results_df(., metadata.table$method[i], metadata.table$dataset_pseudobulk[i],
                       metadata.table$sc_dataset[i])
  
  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table$method[i], metadata.table$dataset_pseudobulk[i],
                       metadata.table$sc_dataset[i], predicted = FALSE)
  
  
  
  result <- left_join(result, true.fractions)
  return(result)
})

data <- bind_rows(tmp)
data <- data[!is.na(data$true_value), ]

# General correlation results
###############################################################################

correlation.results <- data %>%
  group_by(celltype, method, dataset_pseudobulk, dataset_signature) %>%
  dplyr::summarize(corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results.all <- data %>%
  group_by(method, dataset_pseudobulk, dataset_signature) %>%
  dplyr::summarize(celltype = 'all',
                   corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results <- rbind(correlation.results, correlation.results.all)

###############################################################################
# Final figure
###############################################################################

vanderbilt.plot <- correlation.results[correlation.results$celltype!='all', ] %>%
  ggplot(., aes(x=celltype, y=method, fill=corr)) +
  geom_tile()+
  #coord_fixed()+
  geom_text(aes(label = round(corr,2)), size = 2.5)+
  facet_wrap(.~dataset_signature, ncol=1,
             labeller=labeller(
               dataset_signature = c('hao' = 'signature: hao',
                                     'lambrechts' = 'signature: \nlambrechts',
                                     'maynard' = 'signature: \nmaynard')
             ))+
  theme_minimal()+
  scale_fill_gradient2(low = '#ce273f', high='#2d87bb', limits=c(-1,1)) +
  labs(title = 'Vanderbilt', xlab='\ncelltype') +
  theme(axis.title.x = element_text(vjust = -2), legend.position = 'hide') +
  rotate_x_text(angle=60)

ggsave(plot = vanderbilt.plot, '/nfs/proj/omnideconv_benchmarking/omnideconv/benchmark/revision2/visualization/fig6/Fig_6D.pdf', dpi=350, width=2, height=9)
