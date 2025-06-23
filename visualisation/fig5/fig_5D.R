library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
source('./cell_palette.R')
methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c('cpm', rep('counts', 2),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))

# 1: List directories, methods, cell types

unknown.content.deconv.results <- list.files('/vol/omnideconv_results/results_unknown_content', full.names=F, recursive=T)

unknown.content.deconv.results <- unknown.content.deconv.results[grep('replicate', unknown.content.deconv.results)]
unknown.content.deconv.results <- unknown.content.deconv.results[grep('lambrechts', unknown.content.deconv.results)]


metadata.table <- unknown.content.deconv.results %>%
  tibble(path = .,
         method_dataset = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         replicate = gsub('deconvolution__|.rds', '', map_vec(., function(x) strsplit(x, split = '/')[[1]][2])),
         unknown_fraction = map_vec(., function(x) strsplit(x, split = '/')[[1]][3])) %>%
  mutate(method_dataset = gsub("_unknown_content_sim", "", method_dataset),
         replicate = as.numeric(gsub("replicate_", "", replicate)),
         unknown_fraction = as.numeric(gsub("deconvolution_Tumor cells_|.rds", "", unknown_fraction))) %>%
  separate(method_dataset, c("method", "datset"), sep="_")

#metadata.table <- metadata.table[metadata.table$method != 'cibersortx',]


vector.frac <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9)

metadata.table <- metadata.table[metadata.table$unknown_fraction %in% vector.frac, ]



#2: Combine these in a unique dataframe
data <- NULL
for(i in 1:nrow(metadata.table)){
  result <- readRDS(paste0('/vol/omnideconv_results/results_unknown_content/', metadata.table$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  true.fractions <- readRDS(paste0('/vol/omnideconv_results/results_unknown_content/', metadata.table$path[i])) %>%
    .$true_cell_fractions %>%
    as.data.frame() %>%
    .[!(row.names(.) == 'Tumor cells'), ]

  result <- result %>%
    tibble::rownames_to_column(., var='sample') %>%
    gather(., key='celltype', value = 'predicted_value', -'sample') %>%
    mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., method = metadata.table$method[i]) %>%
    mutate(., dataset = metadata.table$datset[i]) %>%
    mutate(., replicate = metadata.table$replicate[i]) %>%
    mutate(., unknown_fraction = metadata.table$unknown_fraction[i])

  true.fractions <- true.fractions %>%
    tibble::rownames_to_column(., var='celltype') %>%
    gather(., key='sample', value = 'true_value', -'celltype') %>%
    mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., method = metadata.table$method[i]) %>%
    mutate(., dataset = metadata.table$datset[i]) %>%
    mutate(., replicate = metadata.table$replicate[i]) %>%
    mutate(., unknown_fraction = metadata.table$unknown_fraction[i])

  result <- left_join(result, true.fractions)
  result$absolute_difference <- abs(result$true_value - result$predicted_value)
  data <- rbind(data, result)
}

# How do we measure how good the deconvolution is? Correlation coefficient

corelation.results <- data %>%
  group_by(celltype, method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value) %>%
  ungroup()

corelation.results.all <- data %>%
  group_by(method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(corr = cor.test(true_value, predicted_value, method = 'pearson')$estimate,
                   pval = cor.test(true_value, predicted_value, method = 'pearson')$p.value) %>%
  ungroup()
corelation.results.all$celltype <- 'all-cells'

rmse.results <- data %>%
  group_by(celltype, method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()

rmse.results.all <- data %>%
  group_by(method, replicate, dataset, unknown_fraction) %>%
  dplyr::summarize(rmse = sqrt(mean((true_value - predicted_value)^2))) %>%
  ungroup()
rmse.results.all$celltype <- 'all-cells'

rmse.results <- rbind(rmse.results, rmse.results.all)
corelation.results <- rbind(corelation.results, corelation.results.all)



data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      median = median(x[[col]], na.rm=TRUE),
      iqr = quantile(x[[col]], probs=c(.25, .75), na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE)
      )
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- plyr::rename(data_sum, c("median" = varname))
  return(data_sum)
}


correlation.results.summary <- data_summary(corelation.results, 'corr', c('celltype', 'method', 'unknown_fraction', 'dataset'))
rmse.results.summary <- data_summary(rmse.results, 'rmse', c('celltype', 'method', 'unknown_fraction', 'dataset'))

# 3: Plot
getPalette = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
colors = getPalette(length(unique(data$celltype)))
names(colors) = unique(data$celltype)

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
                  'all-cells' = "#B3B3B3")

corr.ribbon <- ggplot(correlation.results.summary, aes(x=as.character(unknown_fraction), y=corr, group=celltype))+
  geom_point(aes(color=celltype))+
  geom_line(aes(color=celltype))+
  geom_errorbar(aes(group=unknown_fraction, ymin=`iqr.25%`, ymax=`iqr.75%`), width=.1)+
  geom_ribbon(aes(fill=celltype, ymin=`iqr.25%`, ymax=`iqr.75%`), alpha=.2)+
  facet_wrap(~method, ncol=4)+
  geom_hline(yintercept = 0, linetype = 'dashed')+
  geom_vline(xintercept = '0', linetype = 'dashed') +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  theme_minimal()+
  theme(legend.position = 'hide') +
  #scale_color_brewer(palette = 'Set1')+
  scale_color_manual(values = cell_palette) +
  scale_fill_manual(values = cell_palette) +
  xlab('Unknown cell fraction') +
  ylab('Pearson correlation') +
  labs(title='Correlation trend with increasing tumor fraction')
ggsave(corr.ribbon, "./visualizations/fig_5/fig_5D.pdf", width=13, height = 8)




