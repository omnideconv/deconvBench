library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
library(cowplot)

source('/vol/omnideconv_input/benchmark/pipeline/bin/general_functions/deconvolution_workflow_for_simulation.R')
methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c(rep('counts', 3),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))


##############################################################################
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
##############################################################################

# 1: List directories, methods, cell types
#############################################################################
resolution.deconv.results <- list.files('/vol/omnideconv_results/results_resolution', full.names=F, recursive=T)

resolution.deconv.results <- resolution.deconv.results[grep('wu', resolution.deconv.results)]
resolution.deconv.results <- resolution.deconv.results[grep('deconvolution.rds', resolution.deconv.results)]

metadata.table <- resolution.deconv.results %>%
  tibble(path = .,
         method = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         dataset_level = map_vec(., function(x) strsplit(x, split = '/')[[1]][2]),
         replicate = map_vec(., function(x) strsplit(x, split = '/')[[1]][3])) %>%
  mutate(dataset_level = gsub("_resolution_analysis_sim|_annot", "", dataset_level),
         replicate = gsub('replicate_', '', replicate)) %>%
  separate(dataset_level, c("dataset", "level"), sep="_")


#metadata.table <- metadata.table[metadata.table$method %in% c('scdc', 'cibersortx', 'dwls', 'music', 'bayesprism', 'bisque'), ]
#2: Combine these in a unique dataframe
################################################################################
data <- NULL

table.annotations <- read.table(paste0('/vol/omnideconv_input/omnideconv_data/singleCell/wu', '/cell_type_mappings.csv'), header = T, sep=',')
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
annotation <- read.table(paste0('/vol/omnideconv_input/omnideconv_data/singleCell/wu', '/cell_type_mappings.csv'), header = T, sep=',')
colnames(annotation) <- c('fine', 'normal', 'coarse')

# We read results

for(i in 1:nrow(metadata.table)){
  result <- readRDS(paste0('/vol/omnideconv_results/results_resolution/', metadata.table$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  colnames(result) <- gsub("xxxx", " ", colnames(result))
  colnames(result) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(result))
  colnames(result) <- cleanCelltypesAutogenes(colnames(result))

  true.fractions <- readRDS(paste0('/vol/omnideconv_results/results_resolution/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("xxxx", " ", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))
  colnames(true.fractions) <- cleanCelltypesAutogenes(colnames(true.fractions))


  if(metadata.table$level[i] != 'fine'){
    true.fractions <- reannotate_facs_new(true.fractions, table.annotations, metadata.table$level[i])
  }
  result <- result %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i], metadata.table$replicate[i])

  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i],
                       metadata.table$replicate[i], predicted = FALSE)

  result <- left_join(result, true.fractions)

  data <- rbind(data, result)

}

# We read true fractions

for(i in c(1, 6, 11)){
  true.fractions <- readRDS(paste0('/vol/omnideconv_results/results_resolution/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    t() %>%
    as.data.frame()
  colnames(true.fractions) <- gsub("xxxx", " ", colnames(true.fractions))
  colnames(true.fractions) <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", colnames(true.fractions))
  colnames(true.fractions) <- cleanCelltypesAutogenes(colnames(true.fractions))

  if(metadata.table$level[i] != 'fine'){
    true.fractions <- reannotate_facs_new(true.fractions, table.annotations, metadata.table$level[i])
  }

  #true.fractions <- reannotate_facs_new(true.fractions, annotations, metadata.table$level[i])
  true.fractions <- true.fractions %>%
    process_results_df(., metadata.table$method[i], metadata.table$level[i], 1, predicted = FALSE)

  true.fractions$method <- 'true_values'
  true.fractions$predicted_value <- true.fractions$true_value

  data <- rbind(data, true.fractions)

}

data <- data[data$resolution != 'fine', ]
data$resolution[data$resolution=='normal'] <- 'fine'
data <- data[!is.na(data$true_value), ]
data$resolution <- factor(data$resolution, level=c('coarse', 'fine'))

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

resolution.table <- table.annotations[, -c(1)]
resolution.table$group_coarse <- resolution.table$coarse
resolution.table <- gather(resolution.table, key='res', value='celltype', -group_coarse)
resolution.table$res <- NULL
resolution.table <- distinct(resolution.table)

data.coarsed <- left_join(data, resolution.table)
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


annotation <- table.annotations
table.annotations$fine <- NULL
table.annotations$fine <- table.annotations$normal
table.annotations$normal <- NULL
table.annotations$coarse_celltype <- table.annotations$coarse
annotations <- gather(table.annotations, key='resolution', value='celltype', -coarse_celltype)
annotation$resolution <- NULL

correlation.results <- left_join(correlation.results, annotation)
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
                                'CAFs' = 15,
                                'Myeloid' = 18,
                                'T cells' = 17))

ggsave(plot, filename = './visualisation/supplement/Fig_s6.pdf', width = 10, height = 10)
