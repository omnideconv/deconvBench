source('visualizations/helper_functions.R')
library(tidyverse)
library(data.table)
library(cowplot)

methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c('cpm', rep('counts', 2),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))
method_parameter_df <- data.frame(method=methods, sc_norm=sc_norm, bulk_norm=bulk_norm)
method_parameter_df$method_norm_combi <- paste0(method_parameter_df$method, method_parameter_df$sc_norm, method_parameter_df$bulk_norm)

fractions_df_mouse <- get_fractions('/vol/omnideconv_results/results_mouse/', '0', 'tabula-muris', method_parameter_df)
fractions_df_mouse$org <- 'mm'
fractions_df_hao <- get_all_fractions('/vol/omnideconv_results/results_run/', '0', 'hao-sampled-3', method_parameter_df)
fractions_df_hao$org <- 'hs'
fractions_df <- rbind(fractions_df_hao, fractions_df_mouse)
fractions_df <- subset(fractions_df, method_norm_combi %in% method_parameter_df$method_norm_combi)

performance_df_mouse <- get_performance_metrics('/vol/omnideconv_results/results_mouse/', '0', 'tabula-muris', method_parameter_df)
performance_df_mouse$org <- 'mm'
performance_df_hao <- get_all_performance_metrics('/vol/omnideconv_results/results_run/', '0', 'hao-sampled-3', method_parameter_df)
performance_df_hao$org <- 'hs'
performance_df <- rbind(performance_df_hao, performance_df_mouse)
performance_df <- subset(performance_df, method_norm_combi %in% method_parameter_df$method_norm_combi)

df.cor <- data_summary(performance_df, 'cor', c('cell_type','method','bulk_ds','org'))
df.rmse <- data_summary(performance_df, 'RMSE', c('cell_type','method','bulk_ds','org'))
df.both <- data.table(merge(df.cor, df.rmse, by = c('cell_type','method','bulk_ds','org'), suffixes = c('.cor','.rmse')))

color_palette <- c('B cells'='#999933',
                   'Macrophages'='#CC6677',
                   'mDC'='#882255',
                   'Monocytes'='#AA4499',
                   'NK cells'='#DDCC77',
                   'T cell'='#332288',
                   'T cells CD4 conv'='#117733',
                   'T cells CD8'='#44AA99',
                   'Tregs'='#88CCEE')

df.both$cell_type <- as.factor(df.both$cell_type)
fractions_df$cell_type <- as.factor(fractions_df$celltype)

fractions_df <- fractions_df %>% 
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))
fractions_df <- fractions_df %>%
  mutate(bulk_ds = recode(bulk_ds,
                              "chen" = "Chen",
                              "chen-simulation" = "ChenSim",
                              "finotello" = "Finotello",
                              "finotello-simulation" = "FinotelloSim",
                              "hoek" = "Hoek",
                              "hoek-simulation" = "HoekSim",
                              "petitprez" = "Petitprez",
                              "petitprez-simulation" = "PetitprezSim"))
df.both <- df.both %>%
  mutate(bulk_ds =  recode(bulk_ds,
                           "chen" = "Chen",
                           "chen-simulation" = "ChenSim",
                           "finotello" = "Finotello",
                           "finotello-simulation" = "FinotelloSim",
                           "hoek" = "Hoek",
                           "hoek-simulation" = "HoekSim",
                           "petitprez" = "Petitprez",
                           "petitprez-simulation" = "PetitprezSim"))
df.both <- df.both %>% 
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))

fractions_df <- fractions_df %>%
  mutate(sc_ds = recode(sc_ds,
                          "hao-sampled-3" = "HaoSub",
                          "tabula-muris" = "TM"))


#### Fig 2 ####

fig_2a <- fractions_df %>% subset(bulk_ds %in% c('Hoek')) %>%
  ggplot(., aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(linetype='dashed')+
  geom_point(aes(color=cell_type), size=2.5)+
  facet_wrap( ~ method, ncol=4)+
  xlab('Ground truth cell-type fraction')+
  ylab('Estimated cell-type fraction')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'top',
        strip.background = element_rect(fill = 'white'))+
  scale_color_manual(values = color_palette, drop = T)+
  geom_text(data=df.both %>% subset(bulk_ds %in% c('Hoek') & cell_type == 'all'), 
            aes(label=paste0('R: ',round(cor, 3),'\nRMSE: ', round(RMSE, 3))), x=.05, y=.7, size=3.5, hjust=0)+
  labs(color='Celltype')

fig_2b <- fractions_df %>% subset(bulk_ds %in% c('HoekSim')) %>%
  ggplot(., aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(linetype='dashed')+
  geom_point(aes(color=cell_type), size=2.5)+
  facet_wrap( ~ method, ncol=4)+
  xlab('Ground truth cell-type fraction')+
  ylab('Estimated cell-type fraction')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'top',
        strip.background = element_rect(fill = 'white'))+
  scale_color_manual(values = color_palette, drop = T)+
  geom_text(data=df.both %>% subset(bulk_ds %in% c('HoekSim') & cell_type == 'all'), 
            aes(label=paste0('R: ',round(cor, 3),'\nRMSE: ', round(RMSE, 3))), x=.05, y=.65, size=3.5, hjust=0)+
  labs(color='Celltype')

fig_2c <- df.both %>% subset(cell_type != 'all' & bulk_ds %in% c('Hoek','HoekSim')) %>%
  ggplot(., aes(y=cor, x=RMSE))+
  geom_hline(yintercept = 0, color='#444444', linetype='dashed')+
  geom_hline(yintercept = 0.5, color='#444444', linetype='dashed')+
  geom_hline(yintercept = 1, color='#444444', linetype='dashed')+
  geom_point(aes(color=method), size=2.5)+
  facet_grid(bulk_ds~cell_type)+
  scale_color_brewer(palette = 'Set2')+
  theme_bw()+
  ylab('Pearson Correlation Coefficient')+
  theme(legend.position = 'top',
        strip.background = element_rect(fill = 'white'),
        axis.text.x = element_text(angle = 45, hjust = 1))


ggsave('visualizations_final/fig2/fig_2a.pdf', plot=fig_2a, width = 10, height = 5)
ggsave('visualizations_final/fig2/fig_2b.pdf', plot=fig_2b, width = 10, height = 5)
ggsave('visualizations_final/fig2/fig_2c.pdf', plot=fig_2c, width = 10, height = 5)


#### Fig S2 ####

p.chen.supp <- fractions_df %>% subset(bulk_ds %in% c('Chen','ChenSim')) %>%
  ggplot(., aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(linetype='dashed')+
  geom_point(aes(color=cell_type))+
  facet_grid(bulk_ds ~ method)+
  theme_bw()+
  xlab('Ground truth cell-type fraction')+
  ylab('Estimated cell-type fraction')+
  theme(legend.position = 'top',
        strip.background = element_rect(fill = 'white'))+
  scale_color_manual(values = color_palette, drop = F)+
  geom_text(data=df.both %>% subset(bulk_ds %in% c('Chen','ChenSim') & cell_type == 'all'), 
            aes(label=paste0('R: ',round(cor, 3),'\nRMSE: ', round(RMSE, 3))), x=0, y=.65, size=2.5, hjust=0)

p.peti.supp <- fractions_df %>% subset(bulk_ds %in% c('Petitprez','PetitprezSim')) %>%
  ggplot(., aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(linetype='dashed')+
  geom_point(aes(color=cell_type))+
  facet_grid(bulk_ds ~ method)+
  theme_bw()+
  xlab('Ground truth cell-type fraction')+
  ylab('Estimated cell-type fraction')+
  theme(legend.position = 'top',
        strip.background = element_rect(fill = 'white'))+
  scale_color_manual(values = color_palette, drop = F)+
  geom_text(data=df.both %>% subset(bulk_ds %in% c('Petitprez','PetitprezSim') & cell_type == 'all'), 
            aes(label=paste0('R: ',round(cor, 3),'\nRMSE: ', round(RMSE, 3))), x=0, y=.8, size=2.5, hjust=0)

p.fino.supp <- fractions_df %>% subset(bulk_ds %in% c('Finotello','FinotelloSim')) %>%
  ggplot(., aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(linetype='dashed')+
  geom_point(aes(color=cell_type))+
  facet_grid(bulk_ds ~ method)+
  theme_bw()+
  xlab('Ground truth cell-type fraction')+
  ylab('Estimated cell-type fraction')+
  theme(legend.position = 'top',
        strip.background = element_rect(fill = 'white'))+
  scale_color_manual(values = color_palette, drop = F)+
  geom_text(data=df.both %>% subset(bulk_ds %in% c('Finotello','FinotelloSim') & cell_type == 'all'), 
            aes(label=paste0('R: ',round(cor, 3),'\nRMSE: ', round(RMSE, 3))), x=0, y=.6, size=2.5, hjust=0)+
  labs(color='Celltype')


p.supp <- plot_grid(
  p.fino.supp + theme(legend.position="none"),
  p.peti.supp + theme(legend.position="none"),
  p.chen.supp + theme(legend.position="none"),
  align = 'h',
  labels = c("A", "B", "C"), 
  label_size = 15, 
  hjust = -1,
  nrow = 3, ncol = 1
)

legend <- get_legend(
  p.fino.supp + 
    guides(color = guide_legend(nrow = 1, override.aes = list(size=2))) +
    theme(legend.position = "top")
)

fig_s2 <- plot_grid(
  legend, 
  p.supp, 
  nrow = 2,
  ncol = 1,
  rel_heights = c(.05, .8)
)

ggsave(plot = fig_s2, filename = 'visualizations_final/supplement/fig_s2.png', width = 12, height = 10)


#### Fig S3 ####

fig_s3 <- df.both %>% subset(cell_type != 'all' & bulk_ds %in% c('Finotello','FinotelloSim','Chen','ChenSim','Petitprez','PetitprezSim')) %>%
  ggplot(., aes(y=cor, x=RMSE))+
  geom_hline(yintercept = 0, color='#444444', linetype='dashed')+
  geom_hline(yintercept = 0.5, color='#444444', linetype='dashed')+
  geom_hline(yintercept = 1, color='#444444', linetype='dashed')+
  geom_point(aes(color=method), size=2.5)+
  facet_grid(bulk_ds~cell_type)+
  scale_color_brewer(palette = 'Set2')+
  theme_bw()+
  ylab('Pearson Correlation Coefficient')+
  theme(legend.position = 'top',
        strip.background = element_rect(fill = 'white'),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = fig_s3, filename = 'visualizations_final/supplement/fig_s3.png', width = 12, height = 10)

