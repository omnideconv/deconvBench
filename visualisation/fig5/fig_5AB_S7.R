library(tidyverse)
library(dplyr)
library(ggpubr)
library(gridGraphics)

methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c(rep('counts', 3),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))

# 1: List directories, methods, cell types
missing.cell.types.deconv.results <- list.files('/vol/omnideconv_results/results_missing_cell_types', full.names=F, recursive=T)

metadata.table <- missing.cell.types.deconv.results %>%
  tibble(path = .,
         method = map_vec(., function(x) strsplit(x, split = '/')[[1]][1]),
         metadata = map_vec(., function(x) strsplit(x, split = '/')[[1]][2])) %>%
  mutate(metadata = gsub("no_", "", metadata)) %>%
  mutate(sc_dataset = unlist(lapply(metadata, function(x) stri_split_fixed(x, pattern='_', n=3)[[1]][1])),
         bulk_dataset = unlist(lapply(metadata, function(x) stri_split_fixed(x, pattern='_', n=3)[[1]][2])),
         missing_celltype = unlist(lapply(metadata, function(x) stri_split_fixed(x, pattern='_', n=3)[[1]][3])))

metadata.table$metadata <- NULL
metadata.table <- metadata.table[grepl('deconvolution.rds', metadata.table$path, fixed = TRUE), ]

metadata.table$missing_celltype <- gsub('_', ' ', metadata.table$missing_celltype)


#2: Combine these in a unique dataframe
data <- NULL
for(i in 1:nrow(metadata.table)){
  result <- readRDS(paste0('/vol/omnideconv_results/results_missing_cell_types/', metadata.table$path[i])) %>%
    .$deconv.result %>%
    as.data.frame()

  result <- result %>%
    tibble::rownames_to_column(., var='sample') %>%
    gather(., key='celltype', value = 'predicted_value', -'sample') %>%
    mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., missing_celltype = metadata.table$missing_celltype[i]) %>%
    mutate(., method = metadata.table$method[i]) %>%
    mutate(., bulk_dataset = metadata.table$bulk_dataset[i])
  data <- rbind(data, result)
}

metrics.data <- NULL
for(i in 1:nrow(metadata.table)){
  corr.result <- readRDS(paste0('/vol/omnideconv_results/results_missing_cell_types/', metadata.table$path[i])) %>%
    .$cor_cell_type %>%
    as.data.frame()

  rmse.result <- readRDS(paste0('/vol/omnideconv_results/results_missing_cell_types/', metadata.table$path[i])) %>%
    .$rmse_cell_type %>%
    as.data.frame()

  corr.result <- corr.result %>%
    mutate(., cell_type = gsub("xxxx", " ", cell_type)) %>%
    mutate(., missing_celltype = metadata.table$missing_celltype[i]) %>%
    mutate(., method = metadata.table$method[i]) %>%
    mutate(., bulk_dataset = metadata.table$bulk_dataset[i]) %>%
    select(., -c('pval'))


  rmse.result <- rmse.result %>%
    mutate(., cell_type = gsub("xxxx", " ", cell_type)) %>%
    mutate(., missing_celltype = metadata.table$missing_celltype[i]) %>%
    mutate(., method = metadata.table$method[i]) %>%
    mutate(., bulk_dataset = metadata.table$bulk_dataset[i]) %>%
    select(., -c('zRMSE', 'wRMSE'))
  corr.result <- left_join(corr.result, rmse.result)
  metrics.data <- rbind(metrics.data, corr.result)
}

data <- data %>%
  mutate(method = recode(method,
                         'autogenes' = 'AutoGeneS',
                         'bayesprism' = 'BayesPrism',
                         'bisque' = 'Bisque',
                         'cibersortx' = 'CIBERSORTx',
                         'dwls' = 'DWLS',
                         'music' = 'MuSiC',
                         'scaden' = 'Scaden',
                         'scdc' = 'SCDC'))

metrics.data <- metrics.data %>%
  mutate(method = recode(method,
                         'autogenes' = 'AutoGeneS',
                         'bayesprism' = 'BayesPrism',
                         'bisque' = 'Bisque',
                         'cibersortx' = 'CIBERSORTx',
                         'dwls' = 'DWLS',
                         'music' = 'MuSiC',
                         'scaden' = 'Scaden',
                         'scdc' = 'SCDC'))


data$missing_celltype[data$missing_celltype == 'all-cells'] <- 'None'
metrics.data$missing_celltype[metrics.data$missing_celltype == 'all-cells'] <- 'None'

celltypes_ordered_similarity <- c('Tregs', 'T cells CD4 conv', 'T cells CD8', 'ILC', 'Plasma cells', 'NK cells',
                                  'B cells', 'Monocytes', 'mDC', 'pDC', 'Platelet')


# Difference of RMSEs (final figure)

metrics.all.celltypes <- metrics.data[metrics.data$missing_celltype == 'None', ]
metrics.missing.celltypes <- metrics.data[metrics.data$missing_celltype != 'None', ]

metrics.all.celltypes$RMSE_all_celltypes <- metrics.all.celltypes$RMSE
metrics.all.celltypes[, c('RMSE', 'cor', 'missing_celltype')] <- NULL


metrics.missing.celltypes$RMSE_defective <- metrics.missing.celltypes$RMSE
metrics.missing.celltypes[, c('RMSE', 'cor')] <- NULL

metrics.plot <- left_join(metrics.missing.celltypes,
                          metrics.all.celltypes, multiple = 'all') %>%
  mutate(abs_diff_rmse = RMSE_defective - RMSE_all_celltypes)

finotello.plot <- ggplot(metrics.plot[metrics.plot$bulk_dataset=='finotello', ], aes(x=missing_celltype, y=cell_type, fill=abs_diff_rmse)) +
  geom_tile() +
  coord_equal() +
  geom_text(aes(label = round(abs_diff_rmse,2)), size = 2)+
  facet_wrap(.~method, ncol=3)+
  theme_minimal()+
  scale_fill_gradient2(mid = "white", high = "#8a7495", low="#7f9574", name='RMSE defective - \nRMSE full') +
  rotate_x_text(angle=60) +
  scale_x_discrete(limits=celltypes_ordered_similarity) +
  labs(y="Predicted cell type\n", x="Missing cell type", title='Finotello dataset')

hoek.plot <- ggplot(metrics.plot[metrics.plot$bulk_dataset=='hoek', ], aes(x=missing_celltype, y=cell_type, fill=abs_diff_rmse)) +
  geom_tile() +
  coord_equal() +
  geom_text(aes(label = round(abs_diff_rmse,2)), size = 2)+
  facet_wrap(.~method, ncol=3)+
  theme_minimal()+
  scale_fill_gradient2(mid = "white", high = "#8a7495", low="#7f9574", name='RMSE defective - \nRMSE full') +
  rotate_x_text(angle=60) +
  scale_x_discrete(limits=celltypes_ordered_similarity) +
  labs(y="Predicted cell type\n", x="Missing cell type", title = 'Hoek dataset')


ggsave(finotello.plot, './visualization/fig_5/fig_5A.pdf', dpi = 350, width=12, height=12)
ggsave(hoek.plot, './visualization/supplement/fig_S7.pdf', dpi = 350, width=12, height=12)
