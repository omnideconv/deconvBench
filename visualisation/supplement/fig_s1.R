
library(tidyverse)
library(data.table)
library(cowplot)
library(scattermore)

plot_single_cell_metadata <- function(meta_df, colors, pointsize=0, pixels=c(512,512)){


  count_df <- data.table(meta_df)[,.N,by='cell_type'] %>% arrange(N) %>% mutate(cell_type=factor(cell_type, levels=cell_type))

  plot.composition <- ggplot(count_df, aes(x=cell_type, y=N))+
    geom_segment(aes(x=cell_type, xend=cell_type, y=0, yend=N))+
    geom_point(aes(color=cell_type),size=5) +
    theme_minimal()+
    scale_color_manual(values = colors)+
    coord_flip()+
    xlab('')+
    ylab('count')+
    theme(legend.position = 'none')+
    geom_text(aes(label=N), nudge_y = mean(count_df$N)/3)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = 'None')


  meta_df$cell_type <- factor(meta_df$cell_type, levels = levels(count_df$cell_type))

  plot.genes <- ggplot(meta_df, aes(x=cell_type, y=nFeature_RNA, color=cell_type))+
    geom_violin()+
    coord_flip()+
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 23, fill = "white", color = "black", position = position_dodge(width = 0.75))+
    theme_minimal()+
    xlab('')+
    ylab('Number of expressed genes')+
    scale_color_manual(values = colors)+
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank())+
    theme(axis.text.y = element_text(size = 12),
          panel.grid.major.y = element_blank(),
          legend.position = 'None')


  p.full <- plot_grid(plot.composition + theme(legend.position="none"),
    plot.genes + theme(legend.position="none"),
    align = 'h',
    label_size = 15,
    hjust = -1,
    nrow = 1, ncol = 2, rel_widths = c(0.4, 0.3))

  return(p.full)
}


datasets <- c('TabulaSapiens', 'SchulteSchrepping', 'Lee', 'Heimlich',
              'HaoCleaned', 'DominguezConde', 'BioResourceCollaboration', 'Arunachalam')

base_path <- './singleCellXGeneCensus'
output_path <- './metadata_to_plot'

for (ds in datasets) {
  message("Processing ", ds)

  matrix_file <- file.path(base_path, ds, 'matrix_counts.rds')
  ct_file <- file.path(base_path, ds, 'celltype_annotations.rds')
  batch_file <- file.path(base_path, ds, 'batch.rds')

  full_mat <- readRDS(matrix_file)
  ct <- readRDS(ct_file)
  batch <- readRDS(batch_file)

  meta_df <- data.frame(
    cell_type = ct,
    batch = batch,
    nFeature_RNA = colSums(full_mat > 0)
  )

  rm(full_mat)

  saveRDS(meta_df, file.path(output_path, paste0('meta_', tolower(ds), '.rds')))
}

hao_palette <- c('B cells'='#999933',
                 'Macrophages'='#CC6677',
                 'mDC'='#882255',
                 'Monocytes'='#AA4499',
                 'NK cells'='#DDCC77',
                 'T cell'='#332288',
                 'T cells CD4 conv'='#117733',
                 'T cells CD8'='#44AA99',
                 'Tregs'='#88CCEE',
                 'pDC'='#8d4b00',
                 'Neutrophils'='#fffb6a',
                 'ILC'='#58564f',
                 'Plasma cells'='#969489',
                 'Platelet'='#c2c0b3',
                 'Erythrocytes' = '#332288')



plot_list <- lapply(datasets, function(ds) {
  meta_path <- file.path('/gpfs/gpfs1/scratch/c1041161/TUM_backup/metadata_to_plot', paste0('meta_', tolower(ds), '.rds'))
  meta <- readRDS(meta_path)

  plot_single_cell_metadata(meta, hao_palette)
})

names(plot_list) <- datasets


fig_s1 <- plot_grid(
  plotlist = plot_list[1:8],  # Adjust the range as needed
  labels = c('A','B','C','D','E','F', 'G', 'H'),
  label_size = 15,
  hjust = -1,
  ncol = 2
)


ggsave(plot = fig_s1, filename = '/scratch/c1041161/fig_s1.pdf', height = 15, width = 17)


#### Allen ####
allen_full <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/allen-full/matrix_counts.rds')
allen_ct <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/allen-full/celltype_annotations_normal.rds')
allen_batch <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/allen-full/batch.rds')


meta_df_allen <- data.frame('cell_type' = allen_ct,
                            'batch' = allen_batch,
                            'nFeature_RNA' = colSums(allen_full > 0))

saveRDS(meta_df_allen, '/gpfs/gpfs1/scratch/c1041161/TUM_backup/metadata_to_plot/meta_allen.rds')

#### Maynard ####

maynard_full <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/maynard/matrix_counts.rds')
maynard_ct <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/maynard/celltype_annotations.rds')
maynard_batch <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/maynard/batch.rds')

meta_df_maynard <- data.frame('cell_type' = maynard_ct,
                              'batch' = maynard_batch,
                              'nFeature_RNA' = colSums(maynard_full > 0))


saveRDS(meta_df_allen, '/gpfs/gpfs1/scratch/c1041161/TUM_backup/metadata_to_plot/meta_maynard.rds')

allen_palette <- c('Astro' = '#88CCEE',
                   'Endo' = '#c2c0b3',
                   'Micro' = '#fffb6a',
                   'Inhib' = '#44AA99',
                   'Excit' = '#CC6677',
                   'Oligo' = '#882255',
                   'OPC' = '#AA4499')

rm(maynard_full)
maynard_palette <- c('B cells'='#999933',
                     'Macrophages'='#CC6677',
                     'mDC'='#882255',
                     'Monocytes'='#AA4499',
                     'NK cells'='#DDCC77',
                     'T cells CD4 conv'='#117733',
                     'T cells CD8'='#44AA99',
                     'Tregs'='#88CCEE',
                     'pDC'='#8d4b00',
                     'Neutrophils'='#fffb6a',
                     'Plasma cells'='#969489',
                     "Endothelial cells" = '#c12c2c',
                     "Mast cells" = '#c2c0b3',
                     "Epithelial cell" = '#c1772c',
                     "Tumor cells" = '#0f0b06',
                     "Stromal cells" = '#58564f')

#### Lambrechts ####

lambrechts_full <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/lambrechts/matrix_counts.rds')
lambrechts_ct <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/lambrechts/celltype_annotations.rds')
lambrechts_batch <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/lambrechts/batch.rds')

meta_df_lambrechts <- data.frame('cell_type' = lambrechts_ct,
                                 'batch' = lambrechts_batch,
                                 'nFeature_RNA' = colSums(lambrechts_full > 0))

lambrechts_palette <- c('B cells'='#999933',
                        'Macrophages'='#CC6677',
                        'mDC'='#882255',
                        'Monocytes'='#AA4499',
                        'NK cells'='#DDCC77',
                        'T cells CD4 conv'='#117733',
                        'T cells CD8'='#44AA99',
                        'Tregs'='#88CCEE',
                        'pDC'='#8d4b00',
                        'Neutrophils'='#fffb6a',
                        'Plasma cells'='#969489',
                        "Endothelial cells" = '#c12c2c',
                        "Mast cells" = '#c2c0b3',
                        "Epithelial cell" = '#c1772c',
                        "Tumor cells" = '#0f0b06',
                        "Stromal cells" = '#58564f')
rm(lambrechts_full)
saveRDS(meta_df_lambrechts, '/gpfs/gpfs1/scratch/c1041161/TUM_backup/metadata_to_plot/meta_lambrechts.rds')



#### tabula muris ####

tm_full <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/tabula-muris/matrix_counts.rds')
tm_ct <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/tabula-muris/celltype_annotations.rds')
tm_batch <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/tabula-muris/batch.rds')

meta_df_tm <- data.frame('cell_type' = tm_ct,
                         'batch' = tm_batch,
                         'nFeature_RNA' = colSums(tm_full > 0))

tm_palette <- c('B cells'='#999933',
                'Macrophages'='#CC6677',
                'mDC'='#882255',
                'NK cells'='#DDCC77',
                'T cells CD4 conv'='#117733',
                'T cells CD8'='#44AA99',
                'Tregs'='#88CCEE')

rm(tm_full)
saveRDS(meta_df_tm, '/gpfs/gpfs1/scratch/c1041161/TUM_backup/metadata_to_plot/meta_tm.rds')


#### Wu ####

wu_full <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/wu/matrix_counts.rds')
wu_ct <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/wu/celltype_annotations_coarse.rds')
wu_batch <- readRDS('/gpfs/gpfs1/scratch/c1041161/TUM_backup/singleCell/wu/batch.rds')

meta_df_wu <- data.frame('cell_type' = wu_ct,
                         'batch' = wu_batch,
                         'nFeature_RNA' = colSums(wu_full > 0))
wu_palette <- c('B cells'='#999933',
                'CAFs' = '#333333',
                'Cancer Epithelial' = 'black',
                'Endothelial' = '#c12c2c',
                'Myeloid' = '#AA4499',
                'Normal Epithelial' = '#c1772c',
                'Plasmablasts' = '#969489',
                'PVL' = '#c2c0b3',
                'T-cells' = '#18A749'
)

rm(wu_full)
saveRDS(meta_df_wu, '/gpfs/gpfs1/scratch/c1041161/TUM_backup/metadata_to_plot/wu.rds')


#### full plot ####
plot_list[['Maynard']] <- plot_single_cell_metadata(meta_df_maynard, maynard_palette)
plot_list[['Lambrechts']] <- plot_single_cell_metadata(meta_df_lambrechts, lambrechts_palette)
plot_list[['TabulaMuris']] <- plot_single_cell_metadata(meta_df_tm, tm_palette)
plot_list[['Wu']] <- plot_single_cell_metadata(meta_df_wu, wu_palette)

order = c('HaoCleaned', 'Lambrechts', 'Maynard', 'TabulaMuris', 'Wu',
          'Arunachalam', 'BioResourceCollaboration', 'DominguezConde',
          'Heimlich', 'Lee', 'SchulteSchrepping', 'TabulaSapiens', 'Allen')

plot_list <- plot_list[order]

fig_s1 <- plot_grid(
  plotlist = plot_list,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

labels <- order

n_plots <- 13
ncol <- 2
nrow <- ceiling(n_plots / ncol)


x_coords <- rep(seq(0, 1 - 1/ncol, length.out = ncol), times = nrow)[1:n_plots]

top_y <- 1.00
bottom_y <- 0.15

row_positions <- seq(top_y, bottom_y, length.out = nrow)

y_coords <- rep(row_positions, each = ncol)[1:n_plots]

labels[1] <- 'Hao'

final_plot <- ggdraw(fig_s1) +
  draw_plot_label(labels, x = x_coords + 0.01, y = y_coord, hjust = 0, size = 12)

ggsave(plot = final_plot, filename = './fig_s1_legend.pdf', height = 25, width = 20)
