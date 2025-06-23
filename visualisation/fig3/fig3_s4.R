source("visualisation/helper_functions.R")
library(tidyverse)
library(data.table)
library(cowplot)
#### setup ####

ct_values <- c('5','25','50','100','500','0')
methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c('cpm', rep('counts', 2),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))
method_parameter_df <- data.frame(method=methods, sc_norm=sc_norm, bulk_norm=bulk_norm)
method_parameter_df$method_norm_combi <- paste0(method_parameter_df$method, method_parameter_df$sc_norm, method_parameter_df$bulk_norm)
sc_ds <- 'hao-complete'

method_palette <- palette.colors(palette = "Okabe-Ito")[1:8]
names(method_palette) <- c('AutoGeneS','BayesPrism','Bisque','CIBERSORTx','DWLS','MuSiC','Scaden','SCDC')

#### performance metrics ####
performance_df <- get_all_performance_metrics('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_downsample/', ct_values, 'hao-complete', method_parameter_df)

performance_df <- performance_df %>% 
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))

df.cor <- data_summary(subset(performance_df, method_norm_combi %in% method_parameter_df$method_norm_combi), 'cor', c('ct','cell_type','method','bulk_ds'))
df.cor$ct <- gsub('^0','full', df.cor$ct)
#df.cor$ct <- factor(df.cor$ct, levels = c('5','10','25','50','75','100','300','500','full'))
df.cor$ct <- factor(df.cor$ct, levels = c('5','25','50','100','500','full'))

df.rmse <- data_summary(subset(performance_df, method_norm_combi %in% method_parameter_df$method_norm_combi), 'RMSE', c('ct','cell_type','method','bulk_ds'))
df.rmse$ct <- gsub('^0','full', df.rmse$ct)
#df.rmse$ct <- factor(df.rmse$ct, levels = c('5','10','25','50','75','100','300','500','full'))
df.rmse$ct <- factor(df.rmse$ct, levels = c('5','25','50','100','500','full'))

df.cor <- df.cor %>%
  mutate(bulk_ds =  recode(bulk_ds,
                           "finotello" = "Finotello",
                           "finotello-simulation" = "FinotelloSim",
                           "hoek" = "Hoek",
                           "hoek-simulation" = "HoekSim"))

df.rmse <- df.rmse %>%
  mutate(bulk_ds =  recode(bulk_ds,
                           "finotello" = "Finotello",
                           "finotello-simulation" = "FinotelloSim",
                           "hoek" = "Hoek",
                           "hoek-simulation" = "HoekSim"))


#### Fig 3a ####

df.hoek <- df.rmse %>% subset(bulk_ds %in% c('Hoek','HoekSim'))
df.hoek$bulk_ds <- factor(df.hoek$bulk_ds, levels = c('Hoek','HoekSim'))

tmp <- ggplot(df.hoek, aes(x=ct, y=RMSE, group=method))+
  geom_ribbon(aes(fill=method, ymin=`iqr.25%`, ymax=`iqr.75%`), alpha=.2)+
  geom_line(aes(color=method), alpha=.8)+
  geom_point(data = df.hoek %>% subset(ct == 'full'), size=2, stroke=1)+
  geom_point(aes(color=method), size=1.5, stroke=.6)+
  facet_grid(bulk_ds~cell_type, scales='free_x')+
  theme_bw()+
  ylab('RMSE with ground truth')+
  xlab('number of single cells per cell type')+
  scale_color_manual(values = method_palette)+
  scale_fill_manual(values = method_palette)+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust = .5),
        strip.background = element_rect(fill = 'white'))+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 1, override.aes = list(size=3))) 

ggsave(filename = 'visualizations_final/fig3/fig_3a.pdf', tmp, width = 10, height = 4.5)

fig_3a <- plot_grid(
  tmp,
  align = 'h',
  labels = c("a"), 
  label_size = 15, 
  hjust = -1,
  nrow = 1, ncol = 1
)

#### Fig 3b ####

runtimes_df <- get_runtimes('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_downsample/', 
                            '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/',
                            '/nfs/data/omnideconv_benchmarking_clean/preprocess/',
                            ct_values, sc_ds, methods)
runtimes_df[which(runtimes_df$method %in% c('bayesprism','bisque','music','autogenes') & runtimes_df$process == 'SIGNATURE'),'elapsed'] <- 0

combined_runtimes <- runtimes_df[,sum(elapsed),by=c('method','replicate','ct_fraction','sc_ds','cs_norm','bulk_ds','bulk_norm','n_cells')]
colnames(combined_runtimes)[which(colnames(combined_runtimes) == 'V1')] <- 'elapsed'
combined_runtimes$process <- 'COMBINED'
runtime_df <- rbind(runtimes_df, combined_runtimes, fill=TRUE)

runtime_df$elapsed <- runtime_df$elapsed
runtime_df[which(elapsed == -Inf),'elapsed'] <- 0
runtime_df$n_cells <- as.numeric(runtime_df$n_cells)
#runtime_df$ct_fraction <- log(as.numeric(runtime_df$ct_fraction))

plot_df_runtime <- data_summary(runtime_df, 'elapsed', c('n_cells','process','method','bulk_ds','ct_fraction'))
plot_df_runtime$process <- factor(plot_df_runtime$process, levels = c('SIGNATURE','DECONVOLUTION','COMBINED'))
plot_df_runtime$has_signature <- plot_df_runtime$elapsed > 0

plot_df_runtime$ct_fraction <- gsub('^0','full', plot_df_runtime$ct_fraction)
plot_df_runtime$ct_fraction <- factor(plot_df_runtime$ct_fraction, levels = c('5','10','25','50','75','100','300','500','full'))
plot_df_runtime$elapsed_minutes <- plot_df_runtime$elapsed / 60
plot_df_runtime <- plot_df_runtime %>% subset(bulk_ds %in% c('hoek') & has_signature)

plot_df_runtime <- plot_df_runtime %>% 
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))

fig_3b <- ggplot(plot_df_runtime, aes(x=n_cells, y=elapsed_minutes, group=method))+
    geom_line(aes(color=method))+
    geom_point(data = plot_df_runtime %>% subset(ct_fraction == 'full'), color='black', size=2, stroke=1)+
    geom_point(aes(color=method), size=1.5, stroke=.6)+
    facet_wrap(~process)+
    #geom_ribbon(aes(fill=method, ymin=peak_vmem-sd, ymax=peak_vmem+sd), alpha=.2)+
    theme_bw()+
    ylab('elapsed time (min)')+
    xlab('total number of single cells')+
    scale_color_manual(values = method_palette)+
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust = .5),
          strip.background = element_rect(fill = 'white'))+
    coord_trans(x='log1p',y='log1p')+
    scale_x_continuous(breaks = c(55, 275, 550, 1100, 5000, 153000))

#### Fig 3c ####

df1 <- fread('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/memory_traces/trace_music_dwls_bayesprism_bisque.txt',fill=T, na.strings = c(""))
df2 <- fread('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/memory_traces/trace_autogenes.txt',fill=T, na.strings = c(""))
df3 <- fread('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/memory_traces/trace_scaden.txt',fill=T, na.strings = c(""))
df4 <- fread('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/memory_traces/trace_scdc.txt',fill=T, na.strings = c(""))
df5 <- fread('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/memory_traces/trace_autogenes_music_scdc_fullHao.txt',fill=T, na.strings = c(""))
df6 <- fread('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/memory_traces/trace_bayesprism_fullHao.txt',fill=T, na.strings = c(""))
df <- rbindlist(list(df1,df2,df3,df4,df5,df6))

df <- df[-which(df$task_id == '\t'),]
df <- df[rowSums(is.na(df)) != ncol(df),]

script_lines <- seq(2, nrow(df), by=2)
script_entries <- unlist(df[script_lines,'task_id'])
df <- df[!script_lines]
df$script <- script_entries

df_preprocess <- df %>% subset(process == 'PREPROCESS_SINGLE_CELL')
df_sig <- df %>% subset(process == 'CREATE_SIGNATURE')
df_deconv <- df %>% subset(process == 'DECONVOLUTE')

df_sig$script <- gsub('\'','',df_sig$script)
df_sig$process <- 'SIGNATURE'
df_sig[, c('scriptR','sc_file','sc_anno','batch','sc_ds','sc_norm','bulk_file','bulk_ds','bulk_norm','method','out','bool','rep','ct','cores') := tstrsplit(script, " ")]

df_deconv$script <- gsub('\'','',df_deconv$script)
df_deconv$script <- gsub('\\t',' ',df_deconv$script)
df_deconv$process <- 'DECONVOLUTION'
df_deconv[, c('no','scriptR','sc_file','sc_anno','batch','sc_ds','sc_norm','bulk_file','bulk_ds','bulk_norm','method','out','bool','rep','ct','cores') := tstrsplit(script, " ")]
df_deconv$no <- NULL

df_comb <- rbind(df_deconv, df_sig)

df_comb$process <- factor(df_comb$process, levels = c('SIGNATURE','DECONVOLUTION'))
df_comb$peak_vmem <- unlist(lapply(df_comb$peak_vmem, function(i){
  if(grepl('GB', i)){
    as.numeric(gsub('GB','',i))
  }else{
    as.numeric(gsub('MB','',i)) / 1024
  }
}))


n_cells_df <- rbindlist(lapply(1:nrow(df_comb), function(i){
  ct <- as.numeric(df_comb[i,'ct'])
  rep <- as.numeric(df_comb[i,'rep'])
  path <- paste0('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_downsample/',
                 df_comb[i,'method'],
                 '_hao-complete_',
                 df_comb[i,'sc_norm'],'_',
                 df_comb[i,'bulk_ds'],'_',
                 df_comb[i,'bulk_norm'],'_',
                 'ct', ct, '_rep', rep)
  if(ct == 0){
    sc_data_dir <- paste0('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/', 'hao-complete')
  }else if (as.character(ct-1) %in% c('5','25','50','100','500')){
    sc_data_dir <- paste0('/nfs/data/omnideconv_benchmarking_clean/preprocess/hao-complete_counts_perc',ct-1,'_rep',rep)
  }else{
    return()
  }
  
  cell_anno <- readRDS(paste0(sc_data_dir,'/celltype_annotations.rds'))
  n_cells <- length(cell_anno)
  
  return(list('ct' = as.character(ct),
              'rep' = as.character(rep),
              'n_cells' = n_cells))
}))

df_comb <- unique(join(n_cells_df[,c('ct','n_cells')], df_comb, by='ct'))
df_comb$peak_vmem <- df_comb$peak_vmem
df_comb$rep <- as.numeric(df_comb$rep)
df_comb$n_cells <- as.numeric(df_comb$n_cells)
df_comb[which(method %in% c('bayesprism','bisque','music','autogenes') & process == 'SIGNATURE'),'peak_vmem'] <- 0

plot_df_memory <- data_summary(df_comb, 'peak_vmem', c('ct','n_cells','process','method','bulk_ds'))
plot_df_memory$has_signature <- plot_df_memory$peak_vmem > 0

plot_df_memory$ct <- as.numeric(plot_df_memory$ct)-1
plot_df_memory$ct <- gsub(-1,'full', plot_df_memory$ct)
plot_df_memory$ct <- factor(plot_df_memory$ct, levels = c('5','10','25','50','75','100','300','500','full'))
plot_df_memory <- plot_df_memory %>% subset(bulk_ds %in% c('finotello') &
                                              ct %in% c('5','25','50','100','500','full') &
                                              has_signature)
cx_tmp_df <- data.frame(5,55,'DECONVOLUTION','cibersortx','finotello',NA,NA,NA,FALSE)
colnames(cx_tmp_df) <- colnames(plot_df_memory)
plot_df_memory <- rbind(plot_df_memory, cx_tmp_df)

plot_df_memory <- plot_df_memory %>% 
  mutate(method = recode(method,
                         'autogenes'='AutoGeneS',
                         'bayesprism'='BayesPrism',
                         'bisque'='Bisque',
                         'cibersortx'='CIBERSORTx',
                         'dwls'='DWLS',
                         'music'='MuSiC',
                         'scaden'='Scaden',
                         'scdc'='SCDC'))


fig_3c <- ggplot(plot_df_memory, aes(x=n_cells, y=peak_vmem, group=method))+
  geom_line(aes(color=method))+
  geom_point(data = plot_df_memory %>% subset(ct == 'full'), color='black', size=2, stroke=1)+
  geom_point(aes(color=method), size=1.5, stroke=.6)+
  facet_wrap(~process)+
  #geom_ribbon(aes(fill=method, ymin=peak_vmem-sd, ymax=peak_vmem+sd), alpha=.2)+
  theme_bw()+
  ylab('peak virtual memory (GB)')+
  xlab('total number of single cells')+
  scale_color_manual(values = method_palette)+
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust = .5),
        strip.background = element_rect(fill = 'white'))+
  coord_trans(x='log1p',y='log1p')+
  scale_x_continuous(breaks = c(55, 275, 550, 1100, 5000, 153000))


#### Fig 3 ####

fig_3bc <- plot_grid(
  fig_3b + theme(legend.position="none"),
  fig_3c + theme(legend.position="none"),
  align = 'h',
  labels = c("b", "c"), 
  label_size = 15, 
  hjust = -1,
  nrow = 1, ncol = 2
)

fig_3 <- plot_grid(
  fig_3a, 
  fig_3bc, 
  nrow = 2, 
  rel_heights = c(1, .5)
)

ggsave(filename = 'visualizations_final/fig3/fig_3.pdf', fig_3, width = 12, height = 10)


#### Fig S4a ####

df.fino <- df.rmse %>% subset(bulk_ds %in% c('Finotello','FinotelloSim'))
df.fino$bulk_ds <- factor(df.fino$bulk_ds, levels = c('FinotelloSim','Finotello'))

fig_s4a <- ggplot(df.fino)+
  geom_line(aes(x=ct, y=RMSE, group=method, color=method), alpha=.8)+
  geom_point(aes(x=ct, y=RMSE, group=method, color=method))+
  geom_point(data = df.fino %>% subset(ct == 'full'), mapping = aes(x=ct, y=RMSE, group=method), color='black')+
  geom_ribbon(aes(x=ct, y=RMSE, group=method, fill=method, ymin=`iqr.25%`, ymax=`iqr.75%`), alpha=.2)+
  facet_grid(bulk_ds~cell_type, scales='free_x')+
  theme_bw()+
  ylab('RMSE with ground truth')+
  xlab('number of single cells per cell type')+
  scale_color_manual(values = method_palette)+
  scale_fill_manual(values = method_palette)+
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill = 'white'))


#### Fig S4b ####

fractions_df <- get_fractions('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_downsample/', ct_values, 'hao-complete', method_parameter_df)
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
fractions_df$difference <- fractions_df$fraction.estimate - fractions_df$fraction.true
fractions_df$norm_difference <- fractions_df$difference / (fractions_df$fraction.true + 1e-4)
fractions.summary <- data_summary(fractions_df, 'difference', c('ct','celltype','method', 'bulk_ds'))
fractions.summary$ct <- gsub('^0','full', fractions.summary$ct)
fractions.summary$ct <- factor(fractions.summary$ct, levels = c('5','25','50','100','500','full'))
fractions.summary <- fractions.summary %>%
  mutate(bulk_ds =  recode(bulk_ds,
                           "finotello" = "Finotello",
                           "finotello-simulation" = "FinotelloSim",
                           "hoek" = "Hoek",
                           "hoek-simulation" = "HoekSim"))

fig_s4b <- ggplot(fractions.summary)+
  geom_hline(yintercept = 0, color='#444444', linetype='dashed')+
  #geom_ribbon(aes(x=ct, y=difference, group=method, fill=method, ymin=`iqr.25%`, ymax=`iqr.75%`), alpha=.2)+
  geom_point(aes(x=ct, y=difference, group=method, color=method))+
  geom_point(data = fractions.summary %>% subset(ct == 'full'), mapping = aes(x=ct, y=difference, group=method), color='black')+
  geom_line(aes(x=ct, y=difference, group=method, color=method), alpha=.8)+
  facet_grid(cols = vars(celltype), rows = vars(bulk_ds), scales='free_x')+
  theme_bw()+
  ylab('predicted fraction - true fraction')+
  xlab('number of single cells per cell type')+
  scale_color_manual(values = method_palette)+
  scale_fill_manual(values = method_palette)+
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill = 'white'),
        legend.position = "top")


#### Fig S4 ####

fig_s4 <- plot_grid(
  fig_s4a + theme(legend.position="top"),
  fig_s4b + theme(legend.position="none"),
  labels = c("a","b"), 
  label_size = 15, 
  nrow = 2, 
  rel_heights = c(.7, 1)
)

ggsave(filename = 'visualizations_final/supplement/fig_s4.pdf', plot = fig_s4, width = 12, height = 10)
ggsave(filename = 'visualizations_final/supplement/fig_s4.png', plot = fig_s4, width = 12, height = 10)

