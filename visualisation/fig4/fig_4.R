library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
library(cowplot)
source('./cell_palette.R')
source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')
methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c(rep('counts', 3),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))

# 1: List directories, methods, cell types
#############################################################################
resolution.deconv.results <- list.files('/vol/omnideconv_results/results_resolution', full.names=F, recursive=T)

resolution.deconv.results <- resolution.deconv.results[grep('lambrechts', resolution.deconv.results)]
resolution.deconv.results <- resolution.deconv.results[grep('deconvolution.rds', resolution.deconv.results)]

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
data <- NULL

table.annotations <- read.table(paste0('/vol/omnideconv_input/omnideconv_data/singleCell/lambrechts', '/cell_type_mappings.csv'), header = T, sep=',')
colnames(table.annotations) <- c('fine','normal','coarse')

process_results_df <- function(res, method, resolution, replicate, predicted = TRUE){

  res <- res %>%
    tibble::rownames_to_column(., var='sample') %>%
    gather(., key='celltype', value = 'predicted_value', -'sample') %>%
    #mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., method = method) %>%
    mutate(., resolution = resolution) %>%
    mutate(., replicate = replicate)

  if(!predicted){
    colnames(res)[colnames(res) == 'predicted_value'] <- 'true_value'
  }

  res
}

# We read results

for(i in 1:nrow(metadata.table)){
  result <- readRDS(paste0('/vol/omnideconv_results/results_resolution/', metadata.table$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  colnames(result) <- gsub("xxxx", " ", colnames(result))
  colnames(result) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(result))

  true.fractions <- readRDS(paste0('/vol/omnideconv_results/results_resolution/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))

  result <- result %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i], metadata.table$replicate[i])

  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i],
                       metadata.table$replicate[i], predicted = FALSE)

  result <- left_join(result, true.fractions)

  data <- rbind(data, result)

  #result$absolute_difference <- abs(result$true_value - result$predicted_value)

}

# We read true fractions

for(i in c(1, 6, 11)){
  true.fractions <- readRDS(paste0('/vol/omnideconv_results/results_resolution/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("xxxx", " ", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))

  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i], 1, predicted = FALSE)

  true.fractions$method <- 'true_values'
  true.fractions$predicted_value <- true.fractions$true_value

  data <- rbind(data, true.fractions)

}

# # Lowest res (hierarchical clustering)
# celltypes_ordered_similarity <- c('B cells', 'NK cells', 'T cells NK-like', 'T cells CD8 terminally exhausted',
#                                   'Tregs', 'T cells CD8 naive', 'T cells CD4 non-reg', 'T cells CD8 activated',
#                                   'Monocytes classical', 'Macrophages', 'cDC1', 'cDC2', 'Monocytes non-classical', 'pDCs')
# # All res
# celltypes_levels <- c('B cells', 'NK cells', 'T cells NK-like', 'T cells CD8 terminally exhausted',
#                       'Tregs', 'T cells CD8 naive', 'T cells CD4 non-reg', 'T cells CD8 activated',
#                       'T cells CD4', 'T cells CD8', 'T and NK cells',
#                       'Monocytes classical', 'Macrophages', 'Macrophages-Monocytes',
#                       'cDC1', 'cDC2', 'mDCs', 'Monocytes non-classical', 'Monocytes', 'pDCs', 'DCs')

data <- data%>%
  mutate(celltype_label_plot = recode(celltype,
                                      'Macrophages-Monocytes' = 'Mono/Macro',
                                      'Monocytes' = 'Mono',
                                      'Monocytes classical' = 'Mono class',
                                      'Monocytes non-classical' = 'Mono non-class',
                                      'Macrophages' = 'Macro',
                                      'T cells CD8 activated' = 'T CD8 activated',
                                      'T cells CD8 naive' = 'T CD8 naive',
                                      'T cells CD8 terminally exhausted' = 'T CD8 exhausted',
                                      'T cells NK-like' = 'T NK-like',
                                      'T cells CD8' = 'T CD8',
                                      'T cells CD4' = 'T CD4',
                                      'T cells CD4 non-reg' = 'T CD4 non-reg'))

# We read the resolution table
resolution.table <- read.table('/vol/omnideconv_input/omnideconv_data/singleCell/lambrechts/cell_type_mappings.csv',
                               header=TRUE, sep=',')

resolution.table$group_coarse <- resolution.table$Coarse
resolution.table <- gather(resolution.table, key='res', value='celltype', -group_coarse)
resolution.table$res <- NULL
resolution.table <- distinct(resolution.table)


# Correlation metrics for the results

correlation.results <- data %>%
  group_by(celltype, celltype_label_plot, method, resolution) %>%
  dplyr::summarize(corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results.all <- data %>%
  group_by(method, resolution) %>%
  dplyr::summarize(celltype = 'all',
                   celltype_label_plot='all',
                   corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results <- rbind(correlation.results, correlation.results.all)
correlation.results$resolution <- factor(correlation.results$resolution, levels=c('coarse', 'normal', 'fine'))



# Now we aggregate the results to compute the "coarsed" results

data.coarse.res.all.levels <- left_join(data, resolution.table) %>%
  group_by(sample, replicate, method, group_coarse, resolution) %>%
  dplyr::summarise(
    predicted_value_coarse = sum(predicted_value),
    true_value_coarse = sum(true_value)
  ) %>%
  ungroup()

data.coarse.res.all.levels$celltype <- data.coarse.res.all.levels$group_coarse
data.coarse.res.all.levels$celltype_label_plot <- data.coarse.res.all.levels$group_coarse
data.coarsed <- data.coarse.res.all.levels %>%
  mutate(resolution = recode(resolution,
                             'normal' = 'normal_aggr',
                             'fine' = 'fine_aggr'))
data.coarsed$predicted_value <- data.coarsed$predicted_value_coarse
data.coarsed$true_value <- data.coarsed$true_value_coarse

data.coarsed <- select(data.coarsed, -c('predicted_value_coarse', 'true_value_coarse', 'group_coarse'))
data.coarsed <- data.coarsed[data.coarsed$resolution != 'coarse', ]

data.coarsed$celltype_label_plot <- data.coarsed$celltype
data <- rbind(data.coarsed, data)
data$resolution <- factor(data$resolution, levels=c('fine', 'normal', 'coarse', 'normal_aggr', 'fine_aggr'))

# And the coarsed metrics

correlation.coarse <- data.coarse.res.all.levels %>%
  group_by(celltype_label_plot , group_coarse, method, resolution) %>%
  dplyr::summarize(corr = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$estimate,
                   pval = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value_coarse - predicted_value_coarse)^2))) %>%
  ungroup()

correlation.coarse.all <- data.coarse.res.all.levels %>%
  group_by(method, resolution) %>%
  dplyr::summarize(group_coarse = 'all',
                   celltype_label_plot = 'all',
                   corr = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$estimate,
                   pval = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value_coarse - predicted_value_coarse)^2))) %>%
  ungroup()

correlation.coarse <- rbind(correlation.coarse, correlation.coarse.all)


# We bind with the other results
correlation.coarse$celltype <- correlation.coarse$group_coarse
correlation.coarse$group_coarse <- NULL
correlation.coarse <- correlation.coarse %>%
  mutate(resolution = recode(resolution,
                             'normal' = 'medium_aggr',
                             'fine' = 'fine_aggr'))



correlation.coarse <- correlation.coarse[correlation.coarse$resolution != 'coarse', ]

correlation.results <- rbind(correlation.coarse, correlation.results)
correlation.results$resolution <- factor(correlation.results$resolution, levels=c('fine', 'normal', 'coarse', 'medium_aggr', 'fine_aggr'))
correlation.results <- correlation.results[correlation.results$method != 'true_values', ]

data <- data %>%
  mutate(resolution = recode(resolution,
                             'normal' = 'medium',
                             'normal_aggr' = 'medium_aggr'))

correlation.results <- correlation.results %>%
  mutate(resolution = recode(resolution,
                             'normal' = 'medium'))




table.annotations$coarse_celltype <- table.annotations$coarse

annotations <- gather(table.annotations, key='resolution', value='celltype', -coarse_celltype)

correlation.results <- left_join(correlation.results, annotations[, -c(2)])
correlation.results <- unique(correlation.results)
correlation.results$coarse_celltype[correlation.results$celltype=='all'] <- 'all'

correlation.results <- correlation.results %>%
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))


df <- correlation.results%>% subset(celltype_label_plot != 'all')
df.all <- correlation.results%>% subset(celltype_label_plot == 'all')
plot <- ggplot(df)+
        geom_point(mapping=aes(x=rmse, y=corr, color=celltype_label_plot, shape=coarse_celltype), size=2.5, alpha = .7)+
        geom_point(data = df.all, aes(x=rmse, y=corr), size = 2, color='black', shape=16, alpha = .7)+
        facet_grid(method~resolution)+
        geom_hline(yintercept = 0, linetype='dotted') +
        geom_hline(yintercept = 0.5, linetype='dotted') +
        geom_vline(xintercept = 0.1, linetype='dotted') +
        theme_bw()+
        theme(legend.position = 'bottom',
              strip.background = element_rect(fill = 'white'),
              axis.text = element_text(size = 7)) +
        guides(color=guide_legend(ncol=4),
               shape=guide_legend(ncol=1)) +
        rotate_x_text(angle=60)+
        xlab('RMSE')+ylab('Pearson Correlation')+
        scale_shape_manual(values = c('B cells' = 8,
                                      'DCs' = 15,
                                       'Macrophages-Monocytes' = 18,
                                       'T and NK cells' = 17))

ggsave(plot, filename = './visualizations_final/fig_4/Fig_4.pdf', width = 10, height = 10)

