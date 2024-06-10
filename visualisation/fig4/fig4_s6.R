source("visualisation/helper_functions.R")
library(tidyverse)
library(dplyr)
library(ggpubr)
library(scCustomize)

methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c('cpm', rep('counts', 2),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))
result_path <- '/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_resolution/'


# 1: List directories, methods, cell types
#############################################################################
resolution.deconv.results <- list.files(result_path, full.names=F, recursive=T)

resolution.deconv.results.lambrechts <- resolution.deconv.results[grep('lambrechts', resolution.deconv.results)]
resolution.deconv.results.lambrechts <- resolution.deconv.results.lambrechts[grep('deconvolution.rds', resolution.deconv.results.lambrechts)]

resolution.deconv.results.wu <- resolution.deconv.results[grep('wu', resolution.deconv.results)]
resolution.deconv.results.wu <- resolution.deconv.results.wu[grep('deconvolution.rds', resolution.deconv.results.wu)]


metadata.table.lambrechts <- resolution.deconv.results.lambrechts %>%
  tibble(path = .,
         method = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         dataset_level = map_vec(., function(x) strsplit(x, split = '/')[[1]][2]),
         replicate = map_vec(., function(x) strsplit(x, split = '/')[[1]][3])) %>%
  mutate(dataset_level = gsub("_resolution_analysis_sim|_annot", "", dataset_level),
         replicate = gsub('replicate_', '', replicate)) %>%
  separate(dataset_level, c("dataset", "level"), sep="_")

metadata.table.wu <- resolution.deconv.results.wu %>%
  tibble(path = .,
         method = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         dataset_level = map_vec(., function(x) strsplit(x, split = '/')[[1]][2]),
         replicate = map_vec(., function(x) strsplit(x, split = '/')[[1]][3])) %>%
  mutate(dataset_level = gsub("_resolution_analysis_sim|_annot", "", dataset_level),
         replicate = gsub('replicate_', '', replicate)) %>%
  separate(dataset_level, c("dataset", "level"), sep="_")


###############################################################################
# Lambrechts dataset(fig4)
###############################################################################

data <- NULL

table.annotations <- read.table(paste0('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/lambrechts/', '/cell_type_mappings.csv'), header = T, sep=',')
table.annotations$group_coarse <- table.annotations$Coarse
table.annotations <- gather(table.annotations, key='res', value='celltype', -group_coarse)
table.annotations$res <- NULL
table.annotations <- distinct(table.annotations)


for(i in 1:nrow(metadata.table.lambrechts)){
  result <- readRDS(paste0(result_path, metadata.table.lambrechts$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  colnames(result) <- gsub("xxxx", " ", colnames(result))
  colnames(result) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(result))

  true.fractions <- readRDS(paste0(result_path, metadata.table.lambrechts$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()

  result <- result %>%
    process_resolution_results(., metadata.table.lambrechts$method[i], metadata.table.lambrechts$level[i], metadata.table.lambrechts$replicate[i])

  true.fractions <- true.fractions %>%
    process_resolution_results(., metadata.table.lambrechts$method[i], metadata.table.lambrechts$level[i],
                       metadata.table.lambrechts$replicate[i], predicted = FALSE)

  result <- left_join(result, true.fractions)

  data <- rbind(data, result)

}

data <- data %>%
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))

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

########################

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

data.coarse.res.all.levels <- left_join(data, table.annotations) %>%
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


df <- correlation.results%>% subset(celltype_label_plot != 'all')
df.all <- correlation.results%>% subset(celltype_label_plot == 'all')

custom_pallete <- scCustomize::DiscretePalette_scCustomize(num_colors = 22, palette = 'varibow')

p <- ggplot(df)+
  geom_point(mapping=aes(x=rmse, y=corr, color=celltype, shape=coarse_celltype), size=2.5, alpha = .7)+
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
                                'mDCs' = 15,
                                'pDCs' = 7,
                                'Macrophages-Monocytes' = 18,
                                'T and NK cells' = 17))+
  scale_color_manual(values = custom_pallete)
ggsave(plot = p, filename = 'visualisation/plots/fig4.pdf', width = 9, height = 10)


###############################################################################
# Wu dataset (fig s6)
###############################################################################

data <- NULL

table.annotations <- read.table(paste0('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/', '/cell_type_mappings.csv'), header = T, sep=',')
table.annotations$group_coarse <- table.annotations$Coarse
table.annotations <- gather(table.annotations, key='res', value='celltype', -group_coarse)
table.annotations$res <- NULL
table.annotations <- distinct(table.annotations)



for(i in 1:nrow(metadata.table.wu)){
  result <- readRDS(paste0(result_path, metadata.table.wu$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  colnames(result) <- gsub("xxxx", " ", colnames(result))
  colnames(result) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(result))
  colnames(result) <- cleanCelltypesAutogenes(colnames(result))

  true.fractions <- readRDS(paste0(result_path, metadata.table.wu$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("xxxx", " ", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))
  colnames(true.fractions) <- cleanCelltypesAutogenes(colnames(true.fractions))

  if(metadata.table.wu$level[i] != 'fine'){
    true.fractions <- reannotate_facs(true.fractions, table.annotations, metadata.table.wu$level[i])
  }
  result <- result %>%
    process_resolution_results(., metadata.table.wu$method[i], metadata.table.wu$level[i], metadata.table.wu$replicate[i])

  true.fractions <- true.fractions %>%
    process_resolution_results(., metadata.table.wu$method[i], metadata.table.wu$level[i],
                       metadata.table.wu$replicate[i], predicted = FALSE)

  result <- left_join(result, true.fractions)

  data <- rbind(data, result)

}

# fine resolution will not be used for Wu dataset
data <- data[data$resolution != 'fine', ]
data$resolution[data$resolution=='normal'] <- 'fine'

data <- data %>%
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))

########################
# Correlation metrics for the results

correlation.results <- data %>%
  group_by(celltype, method, resolution) %>%
  dplyr::summarize(corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results.all <- data %>%
  group_by(method, resolution) %>%
  dplyr::summarize(celltype = 'all',
                   corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

correlation.results <- rbind(correlation.results, correlation.results.all)
correlation.results$resolution <- factor(correlation.results$resolution, levels=c('coarse', 'fine'))
correlation.results$celltype <- factor(correlation.results$celltype,
                                       levels = c(unique(correlation.results$celltype)))


# Now we aggregate the results to compute the "coarsed" results

table.annotations <- table.annotations[, -c(1)]
table.annotations$group_coarse <- table.annotations$coarse
table.annotations <- gather(table.annotations, key='res', value='celltype', -group_coarse)
table.annotations$res <- NULL
table.annotations <- distinct(table.annotations)

data.coarsed <- left_join(data, table.annotations)
data.coarse.res.all.levels <- data.coarsed %>%
  group_by(sample, replicate, method, group_coarse, resolution) %>%
  dplyr::summarise(
    predicted_value_coarse = sum(predicted_value),
    true_value_coarse = sum(true_value)
  ) %>%
  ungroup()

data.coarse.res.all.levels$celltype <- data.coarse.res.all.levels$group_coarse


data.coarsed <- data.coarse.res.all.levels %>%
  mutate(resolution = recode(resolution,
                             'fine' = 'fine_aggr'))
data.coarsed$predicted_value <- data.coarsed$predicted_value_coarse
data.coarsed$true_value <- data.coarsed$true_value_coarse

data.coarsed <- select(data.coarsed, -c('predicted_value_coarse', 'true_value_coarse', 'group_coarse'))
data.coarsed <- data.coarsed[data.coarsed$resolution != 'coarse', ]

data <- rbind(data.coarsed, data)

data$resolution <- factor(data$resolution, levels=c('fine', 'coarse', 'fine_aggr'))

# And the coarsed metrics

correlation.coarse <- data.coarse.res.all.levels %>%
  group_by(group_coarse, method, resolution) %>%
  dplyr::summarize(corr = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$estimate,
                   pval = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value_coarse - predicted_value_coarse)^2))) %>%
  ungroup()

correlation.coarse.all <- data.coarse.res.all.levels %>%
  group_by(method, resolution) %>%
  dplyr::summarize(group_coarse = 'all',
                   corr = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$estimate,
                   pval = cor.test(true_value_coarse, predicted_value_coarse, method = 'pearson')$p.value,
                   rmse = sqrt(mean((true_value_coarse - predicted_value_coarse)^2))) %>%
  ungroup()

correlation.coarse <- rbind(correlation.coarse, correlation.coarse.all)


# We bind with the other metrics

correlation.coarse$celltype <- correlation.coarse$group_coarse
correlation.coarse$group_coarse <- NULL
correlation.coarse <- correlation.coarse %>%
  mutate(resolution = recode(resolution,
                             'fine' = 'fine_aggr'))

correlation.coarse <- correlation.coarse[correlation.coarse$resolution != 'coarse', ]

correlation.results <- rbind(correlation.coarse, correlation.results)
correlation.results$resolution <- factor(correlation.results$resolution, levels=c('fine', 'coarse', 'fine_aggr'))
correlation.results <- correlation.results[correlation.results$method != 'true_values', ]


#annotation <- table.annotations
table.annotations$fine <- NULL
table.annotations$fine <- table.annotations$normal
table.annotations$normal <- NULL
table.annotations$coarse_celltype <- table.annotations$coarse
annotations <- gather(table.annotations, key='resolution', value='celltype', -coarse_celltype)
annotations$resolution <- NULL

correlation.results <- left_join(correlation.results, annotations)
correlation.results <- unique(correlation.results)
correlation.results$coarse_celltype[correlation.results$celltype == 'all'] <- 'all'

correlation.results <- correlation.results %>%
  mutate(coarse_celltype = recode(coarse_celltype,
                                  'B-cells' = 'B cells',
                                  'T-cells'='T cells'),
         celltype=recode(celltype,
                         'B-cells' = 'B cells',
                         'T-cells'='T cells'))


df <- correlation.results%>% subset(celltype != 'all')
df.all <- correlation.results%>% subset(celltype == 'all')

custom_pallete <- scCustomize::DiscretePalette_scCustomize(num_colors = 17, palette = 'varibow')

p <- ggplot(df)+
  geom_point(mapping=aes(x=rmse, y=corr, color=celltype, shape=coarse_celltype), size=2.5, alpha = .7)+
  geom_point(data = df.all, aes(x=rmse, y=corr), size = 2.5, color='black', shape=16, alpha = .7)+
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
  scale_color_manual(values = custom_pallete)
ggsave(plot = p, filename = 'visualisation/plots/figs6.pdf', width = 8, height = 10)
