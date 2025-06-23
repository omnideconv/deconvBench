source("visualizations/helper_functions.R")
library(tidyverse)
library(data.table)
library(cowplot)

## get data (subsampling) ##

fractions_df <- get_fractions('/vol/omnideconv_results/results_downsample/', ct_values, 'hao-complete', method_parameter_df)
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
fractions.summary <- data_summary(fractions_df, 'fraction.estimate', c('ct','celltype','sample','fraction.true','method', 'bulk_ds'))
fractions.summary$ct <- gsub('^0','full', fractions.summary$ct)
fractions.summary$ct <- factor(fractions.summary$ct, levels = c('5','25','50','100','500','full'))
performance.all <- performance_df %>% subset(cell_type == 'all')
cor.summary <- data_summary(subset(performance.all, method_norm_combi %in% method_parameter_df$method_norm_combi), 'cor', c('ct','method','bulk_ds'))
rmse.summary <- data_summary(subset(performance.all, method_norm_combi %in% method_parameter_df$method_norm_combi), 'RMSE', c('ct','method','bulk_ds'))
anno_df <- cbind(cor.summary[,c('ct','method','bulk_ds','cor')], rmse.summary[,c('RMSE')])
colnames(anno_df) <- c('ct','method','bulk_ds','cor','RMSE')
anno_df$ct <- gsub('^0','full', anno_df$ct)
anno_df$ct <- factor(anno_df$ct, levels = c('5','25','50','100','500','full'))

fractions.summary <- fractions.summary %>%
  mutate(bulk_ds =  recode(bulk_ds,
                           "finotello" = "Finotello",
                           "finotello-simulation" = "FinotelloSim",
                           "hoek" = "Hoek",
                           "hoek-simulation" = "HoekSim"))
anno_df <- anno_df %>%
  mutate(bulk_ds =  recode(bulk_ds,
                           "finotello" = "Finotello",
                           "finotello-simulation" = "FinotelloSim",
                           "hoek" = "Hoek",
                           "hoek-simulation" = "HoekSim"))

## plot dataset ##

dataset <- 'HoekSim'
hoek_sim <- ggplot(fractions.summary %>% subset(bulk_ds == dataset), aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(alpha=.75, color='#444444', linetype='dashed')+
  geom_point(aes(color=celltype))+
  facet_grid(method~ct)+
  geom_errorbar(aes(group=sample, color=celltype, ymin=(fraction.estimate)-sd, ymax=(fraction.estimate)+sd), width=0, alpha=.8)+
  theme_bw()+
  ylab('estimated cell fractions')+
  xlab('true cell fractions')+
  scale_color_brewer(palette = 'Set2')+
  theme(legend.position = 'top')+
  #coord_cartesian(ylim =c(0, .6), xlim =c(0, .6))+
  geom_text(data=anno_df %>% subset(bulk_ds == dataset), 
            aes(label=paste0('R: ',round(cor, 3),'; RMSE: ', round(RMSE, 3))), x=.25, y=.65, size=2.5)
ggsave(filename = 'visualizations_final/supplement/scatter_hoek_sim.png', hoek_sim, height = 10, width = 10)

dataset <- 'Hoek'
hoek <- ggplot(fractions.summary %>% subset(bulk_ds == dataset), aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(alpha=.75, color='#444444', linetype='dashed')+
  geom_point(aes(color=celltype))+
  facet_grid(method~ct)+
  geom_errorbar(aes(group=sample, color=celltype, ymin=(fraction.estimate)-sd, ymax=(fraction.estimate)+sd), width=0, alpha=.8)+
  theme_bw()+
  ylab('estimated cell fractions')+
  xlab('true cell fractions')+
  scale_color_brewer(palette = 'Set2')+
  theme(legend.position = 'top')+
  #coord_cartesian(ylim =c(0, .6), xlim =c(0, .6))+
  geom_text(data=anno_df %>% subset(bulk_ds == dataset), 
            aes(label=paste0('R: ',round(cor, 3),'; RMSE: ', round(RMSE, 3))), x=.25, y=.65, size=2.5)
ggsave(filename = 'visualizations_final/supplement/scatter_hoek.png', hoek, height = 10, width = 10)

dataset <- 'FinotelloSim'
finotello_sim <- ggplot(fractions.summary %>% subset(bulk_ds == dataset), aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(alpha=.75, color='#444444', linetype='dashed')+
  geom_point(aes(color=celltype))+
  facet_grid(method~ct)+
  geom_errorbar(aes(group=sample, color=celltype, ymin=(fraction.estimate)-sd, ymax=(fraction.estimate)+sd), width=0, alpha=.8)+
  theme_bw()+
  ylab('estimated cell fractions')+
  xlab('true cell fractions')+
  scale_color_brewer(palette = 'Set2')+
  theme(legend.position = 'top')+
  #coord_cartesian(ylim =c(0, .6), xlim =c(0, .6))+
  geom_text(data=anno_df %>% subset(bulk_ds == dataset), 
            aes(label=paste0('R: ',round(cor, 3),'; RMSE: ', round(RMSE, 3))), x=.25, y=.65, size=2.5)
ggsave(filename = 'visualizations_final/supplement/scatter_finotello_sim.png', finotello_sim, height = 10, width = 10)

dataset <- 'Finotello'
finotello <- ggplot(fractions.summary %>% subset(bulk_ds == dataset), aes(x=fraction.true, y=fraction.estimate))+
  geom_abline(alpha=.75, color='#444444', linetype='dashed')+
  geom_point(aes(color=celltype))+
  facet_grid(method~ct)+
  geom_errorbar(aes(group=sample, color=celltype, ymin=(fraction.estimate)-sd, ymax=(fraction.estimate)+sd), width=0, alpha=.8)+
  theme_bw()+
  ylab('estimated cell fractions')+
  xlab('true cell fractions')+
  scale_color_brewer(palette = 'Set2')+
  theme(legend.position = 'top')+
  #coord_cartesian(ylim =c(0, .6), xlim =c(0, .6))+
  geom_text(data=anno_df %>% subset(bulk_ds == dataset), 
            aes(label=paste0('R: ',round(cor, 3),'; RMSE: ', round(RMSE, 3))), x=.25, y=.65, size=2.5)
ggsave(filename = 'visualizations_final/supplement/scatter_finotello.png', finotello, height = 10, width = 10)

