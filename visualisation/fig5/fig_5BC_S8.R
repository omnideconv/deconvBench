library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')

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
                  'Plasma cells'="#FFD92F",
                  'ILC'="#66C2A5",
                  'Platelet'="#B3B3B3")
method_palette <- c('autogenes'="#66C2A5",
                    'bayesprism'="#FC8D62",
                    'bisque'="#8DA0CB",
                    'cibersortx'= "#E78AC3",
                    'dwls'="#A6D854",
                    'music'="#FFD92F",
                    'scaden'="#E5C494",
                    'scdc'= "#B3B3B3")

# 1: List directories, methods, cell types

spillover.deconv.results <- list.files('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_spillover/', full.names=F, recursive=T)

metadata.table <- spillover.deconv.results %>%
  tibble(path = .,
         celltype = gsub('deconvolution_spillover_|.rds', '', map_vec(., function(x) strsplit(x, split = '/')[[1]][2])),
         metadata = map_vec(., function(x) strsplit(x, split = '/')[[1]][1])) %>%
  mutate(metadata = gsub("_spillover_sim", "", metadata)) %>%
  separate(metadata, c("method", "datset", "dataset_rep"), sep="_")

celltypes <- c('B_cells', 'mDC', 'Monocytes', 'NK_cells', 'pDC', 'T_cells_CD4_conv',
               'T_cells_CD8', 'Tregs', 'Macrophages', 'Stromal_cells', 'ILC',
               'Platelet', 'Plasma_cells')

metadata.table <- metadata.table[metadata.table$celltype %in% celltypes, ]
metadata.table$celltype <- gsub('_', ' ', metadata.table$celltype)

metadata.table <- metadata.table[metadata.table$datset %in% c('hao-complete'), ]

metadata.table$dataset_rep <- NULL



#2: Combine these in a unique dataframe
data <- NULL
for(i in 1:nrow(metadata.table)){
  result <- readRDS(paste0('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_spillover/', metadata.table$path[i])) %>%
    .$deconvolution %>%
    as.data.frame()
  result <- result %>%
    tibble::rownames_to_column(., var='sample') %>%
    gather(., key='celltype', value = 'predicted_value', -'sample') %>%
    mutate(., celltype = gsub("xxxx", " ", celltype)) %>%
    mutate(., method = metadata.table$method[i]) %>%
    mutate(., dataset = metadata.table$datset[i]) %>%
    mutate(., true_value = ifelse(celltype == metadata.table$celltype[i], 1, 0), true_celltype = metadata.table$celltype[i])
  data <- rbind(data, result)
}

# figure 5B
###############################################################

data <- arrange(data, dataset, method)

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


celltypes_ordered_similarity <- c('Tregs', 'T cells CD4 conv', 'T cells CD8', 'ILC', 'Plasma cells', 'NK cells',
                                  'B cells', 'Monocytes', 'mDC', 'pDC', 'Platelet')

signals = data[data$dataset == 'hao-complete', ] %>%
  mutate(signal_noise=if_else(celltype == true_celltype, "signal", "noise")) %>%
  group_by(sample, method, true_celltype, signal_noise, dataset) %>%
  dplyr::summarise(estimate=sum(predicted_value)) %>%
  spread(signal_noise, estimate) %>%
  mutate(noise_ratio = noise/(noise+signal)) %>%
  mutate(signal_ratio = signal/(noise+signal)) %>%
  ungroup() %>%
  na.omit()
signal_plot <- ggplot(signals, aes(x=true_celltype, y=signal_ratio, fill=true_celltype)) +
  geom_boxplot(position = position_dodge()) +
  facet_wrap(.~method, labeller = label_wrap_gen(width=10), ncol=1) +
  theme_bw() +
  theme(strip.background=element_rect(colour="black",
                                           fill="white")) +
  labs(y="True cell types fraction\n", x="\nPure celltype") +
  rotate_x_text(angle=60) +
  geom_hline(yintercept=0.5, linetype='dashed') +
  scale_fill_manual(values = cell_palette, name="true celltype") +
  theme(legend.position = "hide") +
  scale_x_discrete(limits=celltypes_ordered_similarity)
ggsave('./visualization/fig_5/fig_5B.pdf', dpi = 350,width = 6, height = 18)


# figure s8
###############################################################

sum.predictions = data[data$dataset == 'hao-complete', ] %>%
  group_by(method, true_celltype, celltype, dataset) %>%
  dplyr::summarise(estimate=sum(predicted_value)) %>%
  select(., -c('dataset'))

sum.predictions$estimate <- sum.predictions$estimate/10
sum.predictions$predicted_celltype <- sum.predictions$celltype

ggplot(sum.predictions, aes(x=true_celltype, y=predicted_celltype, fill=estimate)) +
  geom_tile() +
  geom_text(aes(label = round(estimate,2)), size = 2.5)+
  facet_wrap(.~method)+
  theme_minimal()+
  theme(strip.text = element_text(face = "bold")) +
  rotate_x_text(angle = 60) +
  scale_fill_gradient(low='white',high = 'purple', name='Mean estimate') +
  xlab('True celltype') +
  ylab('Predicted celltype') +
  scale_x_discrete(limits=celltypes_ordered_similarity) +
  scale_y_discrete(limits=celltypes_ordered_similarity)

ggsave('./visualization/supplement/fig_S8.pdf', dpi = 350, width = 11, height = 9)

# figure 5C
###############################################################

### code adopted from immunedeconv benchmarking / spillover analysis

colors = cell_palette
names(colors) = names(cell_palette)



resultDfIn = data[data$dataset=='hao-complete', ]
overviewTable = data[data$dataset=='hao-complete', ]
par(mar=rep(2, 4))
circos.par(cell.padding = rep(2, 4))
pdf('./visualizations/fig_5/fig_5C.pdf', width = 22, height = 4)
layout(t(matrix(seq(1, length(unique(overviewTable$dataset)) * length(unique(overviewTable$method))),
                length(unique(overviewTable$method)),
                1)))

for(dataset in unique(data$dataset)){

  for(method in unique(data$method)){
    resultDf <- data
    migration = resultDf %>%
        filter(method == !!method, dataset == !!dataset) %>%
        group_by(method, true_celltype, celltype) %>%
        dplyr::summarise(estimate = mean(predicted_value)) %>%
        ungroup()
    migration_mat = migration %>%
        select(-method) %>%
        spread(celltype, estimate) %>%
        as.data.frame() %>%
        tibble::column_to_rownames("true_celltype") %>%
        as.matrix()
    noise_ratio = migration %>%
        mutate(type = ifelse(celltype == true_celltype, "signal", "noise")) %>%
        group_by(method, type) %>%
        dplyr::summarise(estimate = sum(estimate)) %>%
        spread(type, estimate) %>%
        mutate(noise_ratio = noise/(signal+noise), signal_ratio = signal/(signal+noise)) %>%
        ungroup()
    par(cex = 0.9, mar = c(1, 1, 1, 1))
    chordDiagram(migration_mat, directional = TRUE, transparency = .5,
                   grid.col = colors,
                   annotationTrack =  c("name", "grid")
                   #annotationTrack = c("grid", "axis"), annotationTrackHeight = mm_h(5)
      )
      circos.clear()
      for(si in get.all.sector.index()) {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
        circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1,
                    facing = "bending.inside", niceFacing = TRUE, col = "white")
      }
      text(0, 0, method, cex = 1.8)
      text(0, -0.2, as.character(round(filter(noise_ratio, method == !!method) %>% pull(noise_ratio), 2)), cex=1.8)

    }
  }
dev.off()
