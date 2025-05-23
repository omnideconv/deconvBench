library(Seurat)
library(tidyverse)
library(data.table)
library(cowplot)
library(scattermore)

plot_single_cell <- function(obj, colors, pointsize=0, pixels=c(512,512)){
  meta_df <- obj@meta.data
  meta_df$umap1 <- obj@reductions$umap@cell.embeddings[,1]
  meta_df$umap2 <- obj@reductions$umap@cell.embeddings[,2]
  
  plot.umap <- ggplot(meta_df, aes(x=umap1, y=umap2, color=cell_type))+
    geom_scattermore(aes(color=cell_type), pointsize = pointsize, pixels = pixels)+
    theme_minimal()+
    xlab('UMAP1')+ylab("UMAP2")+
    scale_color_manual(values = colors)+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5)
    )+
    guides(color = guide_legend(title = '',override.aes = list(size = 5)))
  
  
  count_df <- data.table(meta_df)[,.N,by='cell_type'] %>% arrange(N) %>% mutate(cell_type=factor(cell_type, levels=cell_type))
  
  plot.composition <- ggplot(count_df, aes(x=cell_type, y=N))+
    geom_segment(aes(x=cell_type, xend=cell_type, y=0, yend=N))+
    geom_point(aes(color=cell_type),size=5) + 
    theme_minimal()+
    scale_color_manual(values = colors)+
    coord_flip()+
    xlab('')+ylab('count')+
    theme(legend.position = 'none')+
    geom_text(aes(label=N), nudge_y = mean(count_df$N)/5)+
    theme(axis.text.y = element_text(size = 12),
          panel.grid.major.y = element_blank())
  
  
  meta_df$cell_type <- factor(meta_df$cell_type, levels = levels(count_df$cell_type))
  
  plot.genes <- ggplot(meta_df, aes(x=cell_type, y=nFeature_RNA, color=cell_type))+
    geom_violin()+
    coord_flip()+
    stat_summary(fun = "mean", geom = "point", size = 3, shape = 23, fill = "white", color = "black", position = position_dodge(width = 0.75))+
    theme_minimal()+
    xlab('')+ylab('Number of expressed genes')+
    scale_color_manual(values = colors)+
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank())+
    theme(legend.position = 'none')
  
  
  p.full <- plot_grid(
    plot.umap + theme(legend.position="none"),
    plot.composition + theme(legend.position="none"),
    plot.genes + theme(legend.position="none"),
    align = 'h',
    label_size = 15, 
    hjust = -1,
    nrow = 1, ncol = 3, rel_widths = c(.3, 0.4, .3)
  )
  
  return(p.full)
}

#### Hao ####
hao_full <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/matrix_counts.rds')
hao_ct <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/celltype_annotations.rds')
hao_batch <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/batch.rds')

hao_subset_3_cells <- colnames(readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-sampled-3/matrix_counts.rds'))

meta_df <- data.frame('cell_type' = hao_ct,
                      'batch' = hao_batch,
                      'barcode' = colnames(hao_full),
                      row.names = colnames(hao_full))
meta_df$subset <- 'full dataset'
meta_df$subset <- ifelse(meta_df$barcode %in% hao_subset_3_cells, 'subset', meta_df$subset)

hao <- CreateSeuratObject(hao_full, project='hao', meta.data = meta_df)
hao <- SCTransform(hao) %>%
  RunPCA(npcs = 20, features = VariableFeatures(hao)) %>%
  RunUMAP(reduction = "pca", dims = 1:20) 
saveRDS(hao, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/seurat_obj.rds')
hao <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/seurat_obj.rds')

#hao_subset <- subset(hao, subset=='subset')

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
                  'Platelet'='#c2c0b3')

#### Maynard ####
maynard_full <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/maynard/matrix_counts.rds')
maynard_ct <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/maynard/celltype_annotations.rds')
maynard_batch <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/maynard/batch.rds')

meta_df <- data.frame('cell_type' = maynard_ct,
                      'batch' = maynard_batch,
                      'barcode' = colnames(maynard_full),
                      row.names = colnames(maynard_full))
maynard <- CreateSeuratObject(maynard_full, project='maynard', meta.data = meta_df)
maynard <- SCTransform(maynard) %>%
  RunPCA(npcs = 20, features = VariableFeatures(maynard)) %>%
  RunUMAP(reduction = "pca", dims = 1:20) 
saveRDS(maynard, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/maynard/seurat_obj.rds')
maynard <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/maynard/seurat_obj.rds')

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

lambrechts_full <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/lambrechts/matrix_counts.rds')
lambrechts_ct <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/lambrechts/celltype_annotations.rds')
lambrechts_batch <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/lambrechts/batch.rds')

meta_df <- data.frame('cell_type' = lambrechts_ct,
                      'batch' = lambrechts_batch,
                      'barcode' = colnames(lambrechts_full),
                      row.names = colnames(lambrechts_full))
lambrechts <- CreateSeuratObject(lambrechts_full, project='lambrechts', meta.data = meta_df)
lambrechts <- SCTransform(lambrechts) %>%
  RunPCA(npcs = 20, features = VariableFeatures(lambrechts)) %>%
  RunUMAP(reduction = "pca", dims = 1:20) 
saveRDS(lambrechts, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/lambrechts/seurat_obj.rds')
lambrechts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/lambrechts/seurat_obj.rds')

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

#### tabula muris ####

tm_full <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/matrix_counts.rds')
tm_ct <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/celltype_annotations.rds')
tm_batch <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/batch.rds')

meta_df <- data.frame('cell_type' = tm_ct,
                      'batch' = tm_batch,
                      'barcode' = colnames(tm_full),
                      row.names = colnames(tm_full))
tm <- CreateSeuratObject(tm_full, project='tabula-muris', meta.data = meta_df)
tm <- SCTransform(tm) %>%
  RunPCA(npcs = 20, features = VariableFeatures(tm)) %>%
  RunUMAP(reduction = "pca", dims = 1:20) 
saveRDS(tm, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/seurat_obj.rds')
tm <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/seurat_obj.rds')

tm_palette <- c('B cells'='#999933',
                'Macrophages'='#CC6677',
                'mDC'='#882255',
                'NK cells'='#DDCC77',
                'T cells CD4 conv'='#117733',
                'T cells CD8'='#44AA99',
                'Tregs'='#88CCEE'
                )

#### Wu ####

wu_full <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/matrix_counts.rds')
wu_ct <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/celltype_annotations_coarse.rds')
wu_batch <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/batch.rds')

meta_df <- data.frame('cell_type' = wu_ct,
                      'batch' = wu_batch,
                      'barcode' = colnames(wu_full),
                      row.names = colnames(wu_full))
wu <- CreateSeuratObject(wu_full, project='Wu', meta.data = meta_df)
wu <- SCTransform(wu) %>%
  RunPCA(npcs = 20, features = VariableFeatures(wu)) %>%
  RunUMAP(reduction = "pca", dims = 1:20) 
saveRDS(wu, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/seurat_obj.rds')
wu <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/wu/seurat_obj.rds')

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


#### full plot ####

plots.hao <- plot_single_cell(hao, hao_palette, pointsize = 1, pixels = c(1024, 1024))
plots.maynard <- plot_single_cell(maynard, maynard_palette, pointsize = 1.5)
plots.lambrechts <- plot_single_cell(lambrechts, lambrechts_palette, pointsize = 1.5)
plots.tm <- plot_single_cell(tm, tm_palette, pointsize = 2)
plots.wu <- plot_single_cell(wu, wu_palette, pointsize = 1.5)

fig_s1 <- plot_grid(
  plots.hao,
  plots.maynard,
  plots.lambrechts,
  plots.tm,
  plots.wu,
  align = 'h', 
  labels = c('A','B','C','D','E'),
  label_size = 15, 
  hjust = -1,
  nrow = 4, ncol = 1
)

ggsave(plot = fig_s1, filename = 'visualisation/plots/fig_s1/fig_s1.pdf', height = 21, width = 14)
ggsave(plot = plots.wu, filename = 'visualisation/plots/fig_s1/fig_s1_wu.pdf', height = 4.5, width = 14)
