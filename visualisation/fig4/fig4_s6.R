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
    process_results_df(., metadata.table.lambrechts$method[i], metadata.table.lambrechts$level[i], metadata.table.lambrechts$replicate[i])
  
  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table.lambrechts$method[i], metadata.table.lambrechts$level[i], 
                       metadata.table.lambrechts$replicate[i], predicted = FALSE)
  
  result <- left_join(result, true.fractions)
  
  data <- rbind(data, result)
  
}

########################

# calc aggregated values
data.coarsed <- left_join(data, table.annotations)
data.coarse.res.all.levels <- data.coarsed %>%
  group_by(sample, replicate, method, group_coarse, resolution) %>%
  dplyr::summarise(
    predicted_value_coarse = sum(predicted_value),
    true_value_coarse = sum(true_value)
  ) %>%
  ungroup()

data.coarse.res.all.levels$celltype <- data.coarse.res.all.levels$group_coarse
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

# calc resolution-specific values
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


corr.coarse.bind <- correlation.coarse
corr.coarse.bind$celltype <- corr.coarse.bind$group_coarse
corr.coarse.bind$group_coarse <- NULL
corr.coarse.bind <- corr.coarse.bind %>%
  mutate(resolution = recode(resolution,
                             'normal' = 'medium_aggr',
                             'fine' = 'fine_aggr'))
corr.coarse.bind <- corr.coarse.bind[corr.coarse.bind$resolution != 'coarse', ]
corr.coarse.bind <- rbind(corr.coarse.bind, correlation.results)
corr.coarse.bind$resolution <- factor(corr.coarse.bind$resolution, levels=c('fine', 'normal', 'coarse', 'medium_aggr', 'fine_aggr'))
corr.coarse.bind <- corr.coarse.bind[corr.coarse.bind$method != 'true_values', ]

table.annotations <- read.table(paste0('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/lambrechts/', '/cell_type_mappings.csv'), header = T, sep=',')
colnames(table.annotations) <- c('fine','normal','coarse')
table.annotations$coarse_celltype <- table.annotations$coarse
annotation.2 <- gather(table.annotations, key='resolution', value='celltype', -coarse_celltype)
corr.coarse.bind.2 <- left_join(corr.coarse.bind, annotation.2[, -c(2)])
corr.coarse.bind.2 <- unique(corr.coarse.bind.2)


corr.coarse.bind.2 <- corr.coarse.bind.2 %>%
  mutate(celltype = recode(celltype,
                           'Macrophages-Monocytes' = 'Mono/Macro',
                           'Monocytes' = 'Mono',
                           'Monocytes classical' = 'Mono cl',
                           'Monocytes non-classical' = 'Mono ncl',
                           'Macrophages' = 'Macro',
                           'T cells CD8 activated' = 'CD8 act',
                           'T cells CD8 naive' = 'CD8 naiv',
                           'T cells CD8 terminally exhausted' = 'CD8 exh',
                           'T cells NK-like' = 'NK T',
                           'NK cells' = 'NK',
                           'T cells CD8' = 'CD8 T',
                           'T cells CD4' = 'CD4 T',
                           'T cells CD4 non-reg' = 'CD4 non-reg'))

custom_pallete <- scCustomize::DiscretePalette_scCustomize(num_colors = 21, palette = 'varibow')

df <- corr.coarse.bind.2%>% subset(celltype != 'all')
df2 <- corr.coarse.bind.2%>% subset(celltype == 'all')
p <- ggplot(df)+
  geom_point(mapping=aes(x=rmse, y=corr, color=celltype, shape=coarse_celltype), size=2.5, alpha = .7)+
  geom_point(data = df2, aes(x=rmse, y=corr), size = 2, color='black', shape=16, alpha = .7)+
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
                                'T and NK cells' = 17))+
  scale_color_manual(values = custom_pallete)
ggsave(plot = p, filename = 'visualisation/plots/fig4.pdf', width = 9, height = 10)


###############################################################################
# Wu dataset (fig s6)
# ugly code repetition, but works :`)
###############################################################################

data <- NULL

table.annotations <- read.table(paste0('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/', '/cell_type_mappings.csv'), header = T, sep=',')
table.annotations$group_coarse <- table.annotations$Coarse
table.annotations <- gather(table.annotations, key='res', value='celltype', -group_coarse)
table.annotations$res <- NULL
table.annotations <- distinct(table.annotations)

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

cleanCelltypesAutogenes <- function(celltype){
  celltype <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "\\+", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9f", "\\(", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9g", ")", celltype)
  return(gsub("xxxx", " ", celltype))
}

reannotate_facs_new <- function(facs.table, annotation, new.annotation.level){
  
  facs.table <- as.data.frame(facs.table)
  annotation <- annotation[which(annotation$fine %in% colnames(facs.table)), ]
  cell_types <- unique(annotation[[new.annotation.level]])
  for(c in cell_types){
    # These are the fine cell types
    cur.cell.types <- annotation[which(annotation[[new.annotation.level]] == c), 1]
    
    if(length(cur.cell.types) > 1){
      facs.table[c] <- rowSums(facs.table[, c(cur.cell.types)])
      facs.table[, c(cur.cell.types)] <- NULL
    } else if(length(cur.cell.types) == 1) {
      if(c != cur.cell.types){
        facs.table[c] <- facs.table[cur.cell.types]
        facs.table[cur.cell.types] <- NULL
      }
    }
  }
  facs.table
}

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
    true.fractions <- reannotate_facs_new(true.fractions, table.annotations, metadata.table.wu$level[i])
  }
  result <- result %>%
    process_results_df(., metadata.table.wu$method[i], metadata.table.wu$level[i], metadata.table.wu$replicate[i])
  
  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table.wu$method[i], metadata.table.wu$level[i], 
                       metadata.table.wu$replicate[i], predicted = FALSE)
  
  result <- left_join(result, true.fractions)
  
  data <- rbind(data, result)

}

# fine resolution will not be used for Wu dataset
data <- data[data$resolution != 'fine', ]
data$resolution[data$resolution=='normal'] <- 'fine'


########################

# calc aggregated values
data.coarsed <- left_join(data, table.annotations)
data.coarse.res.all.levels <- data.coarsed %>%
  group_by(sample, replicate, method, group_coarse, resolution) %>%
  dplyr::summarise(
    predicted_value_coarse = sum(predicted_value),
    true_value_coarse = sum(true_value)
  ) %>%
  ungroup()

data.coarse.res.all.levels$celltype <- data.coarse.res.all.levels$group_coarse
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

# calc resolution-specific values
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


corr.coarse.bind <- correlation.coarse
corr.coarse.bind$celltype <- corr.coarse.bind$group_coarse
corr.coarse.bind$group_coarse <- NULL
corr.coarse.bind <- corr.coarse.bind %>%
  mutate(resolution = recode(resolution,
                             'fine' = 'fine_aggr'))

corr.coarse.bind <- corr.coarse.bind[corr.coarse.bind$resolution != 'coarse', ]
corr.coarse.bind <- rbind(corr.coarse.bind, correlation.results)
corr.coarse.bind$resolution <- factor(corr.coarse.bind$resolution, levels=c('fine', 'coarse', 'fine_aggr'))
corr.coarse.bind <- corr.coarse.bind[corr.coarse.bind$method != 'true_values', ]

table.annotations <- read.table(paste0('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/', '/cell_type_mappings.csv'), header = T, sep=',')
colnames(table.annotations) <- c('fine','normal','coarse')
table.annotations$fine <- NULL
table.annotations$fine <- table.annotations$normal
table.annotations$normal <- NULL
table.annotations$coarse_celltype <- table.annotations$coarse
annotation.2 <- gather(table.annotations, key='resolution', value='celltype', -coarse_celltype)
annotation.2$resolution <- NULL

corr.coarse.bind.2 <- left_join(corr.coarse.bind, annotation.2)
corr.coarse.bind.2 <- unique(corr.coarse.bind.2)
corr.coarse.bind.2$coarse_celltype[corr.coarse.bind.2$celltype == 'all'] <- 'all'

corr.coarse.bind.2 <- corr.coarse.bind.2 %>%
  mutate(coarse_celltype = recode(coarse_celltype,
                                  'B-cells' = 'B cells',
                                  'T-cells'='T cells'), 
         celltype=recode(celltype,
                         'B-cells' = 'B cells',
                         'T-cells'='T cells'))

custom_pallete <- scCustomize::DiscretePalette_scCustomize(num_colors = 17, palette = 'varibow')

df <- corr.coarse.bind.2%>% subset(celltype != 'all')
df2 <- corr.coarse.bind.2%>% subset(celltype == 'all')
p <- ggplot(df)+
  geom_point(mapping=aes(x=rmse, y=corr, color=celltype, shape=coarse_celltype), size=2.5, alpha = .7)+
  geom_point(data = df2, aes(x=rmse, y=corr), size = 2.5, color='black', shape=16, alpha = .7)+
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
