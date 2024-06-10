library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
library(cowplot)

source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')
methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c(rep('counts', 3),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))

# 1: List directories, methods, cell types
#############################################################################
impact.technology.results <- list.files('/vol/omnideconv_results/results_impact_technology', full.names=F, recursive=T)

#impact.technology.results <- impact.technology.results[grep('lambrechts', impact.technology.results)]
impact.technology.results <- impact.technology.results[grep('deconvolution.rds', impact.technology.results)]

metadata.table <- impact.technology.results %>%
  tibble(path = .,
         method_pseudobulk = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         pseudobulk_signature = map_vec(., function(x) strsplit(x, split = '/')[[1]][2]),
         replicate = map_vec(., function(x) strsplit(x, split = '/')[[1]][3])) %>%
  separate(method_pseudobulk, c("method", "dataset_pseudobulk"), sep="_") %>%
  mutate(dataset_signature = map_vec(pseudobulk_signature, function(x) strsplit(x, split = '_')[[1]][3]),
         replicate = gsub('replicate_', '', replicate))

metadata.table$pseudobulk_signature <- NULL


#2: Combine these in a unique dataframe
################################################################################
data <- NULL

process_results_df <- function(res, method, pseudobulk_ds, signature_ds, replicate, predicted = TRUE){

  res <- res %>%
    tibble::rownames_to_column(., var='sample') %>%
    gather(., key='celltype', value = 'predicted_value', -'sample') %>%
    mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., method = method) %>%
    mutate(., dataset_pseudobulk = pseudobulk_ds) %>%
    mutate(., dataset_signature = signature_ds) %>%
    mutate(., replicate = replicate)

  if(!predicted){
    colnames(res)[colnames(res) == 'predicted_value'] <- 'true_value'
  }

  res
}

for(i in 1:nrow(metadata.table)){
  result <- readRDS(paste0('/vol/omnideconv_results/results_impact_technology/', metadata.table$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  colnames(result) <- gsub("xxxx", " ", colnames(result))
  colnames(result) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(result))

  true.fractions <- readRDS(paste0('/vol/omnideconv_results/results_impact_technology/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("xxxx", " ", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))

  result <- result %>%
    process_results_df(., metadata.table$method[i], metadata.table$dataset_pseudobulk[i],
                       metadata.table$dataset_signature[i], metadata.table$replicate[i])

  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table$method[i], metadata.table$dataset_pseudobulk[i],
                       metadata.table$dataset_signature[i], metadata.table$replicate[i], predicted = FALSE)

  result <- left_join(result, true.fractions)

  data <- rbind(data, result)

  #result$absolute_difference <- abs(result$true_value - result$predicted_value)

}



# General correlation results
###############################################################################
data <- data[data$dataset_pseudobulk %in% c('lambrechts', 'maynard'), ]
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

lambrechts.plot <- ggplot(correlation.results[correlation.results$dataset_pseudobulk == 'lambrechts', ],
                          aes(x=celltype, y=method, fill=corr)) +
  geom_tile()+
  coord_fixed()+
  geom_text(aes(label = round(corr,2)), size = 2.5)+
  facet_wrap(.~dataset_signature, ncol=3,
             labeller=labeller(
               dataset_signature = c('hao-complete' = 'signature: hao',
                                     'lambrechts' = 'signature: lambrechts',
                                     'maynard' = 'signature: maynard')
             ))+
  theme_minimal()+
  scale_fill_gradient2(low = '#ce273f', high='#2d87bb', limits=c(-1,1)) +
  labs(title = 'Lambrechts', xlab='\ncelltype') +
  theme(axis.title.x = element_text(vjust = -2), legend.position = 'hide') +
  rotate_x_text(angle=60)

ggsave(lambrechts.plot, './visualisation/fig_6/fig_6A.pdf', dpi = 350, width=10, height=5)

maynard.plot <- ggplot(correlation.results[correlation.results$dataset_pseudobulk == 'maynard', ],
                       aes(x=celltype, y=method, fill=corr)) +
  geom_tile()+
  coord_fixed()+
  geom_text(aes(label = round(corr,2)), size = 2.5)+
  facet_wrap(.~dataset_signature, ncol=3,
             labeller=labeller(
               dataset_signature = c('hao-complete' = 'signature: hao',
                                     'lambrechts' = 'signature: lambrechts',
                                     'maynard' = 'signature: maynard')
             ))+
  theme_minimal()+
  scale_fill_gradient2(low = '#ce273f', high='#2d87bb', limits=c(-1,1)) +
  labs(title = 'Maynard', xlab='\n\ncelltype\n') +
  theme(axis.title.x = element_text(vjust = -1), legend.position = 'bottom') +
  rotate_x_text(angle=60)


ggarrange(lambrechts.plot,
          NULL,
          maynard.plot,
          nrow = 3, heights = c(7, 0.2, 7))

ggsave(maynard.plot, './visualisation/fig_6/fig_6B.pdf', dpi = 350, width=10, height=5)
