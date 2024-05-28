source("visualisation/helper_functions.R")
library(tidyverse)
library(data.table)
library(cowplot)
#### setup ####

methods <- c('autogenes','bayesprism','bisque','cibersortx','dwls','music','scaden','scdc')
sc_norm <- c('cpm', rep('counts', 2),'cpm',rep('counts',4))
bulk_norm <- c('tpm', rep('counts', 2), rep('tpm',5))
method_parameter_df <- data.frame(method=methods, sc_norm=sc_norm, bulk_norm=bulk_norm)
method_parameter_df$method_norm_combi <- paste0(method_parameter_df$method, method_parameter_df$sc_norm, method_parameter_df$bulk_norm)

performance_df_mouse <- get_performance_metrics('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_mouse/', '0', 'tabula-muris', method_parameter_df)
performance_df_mouse$org <- 'mm'
performance_df_hao <- get_all_performance_metrics('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_main/', '0', 'hao-sampled-3', method_parameter_df)
performance_df_hao$org <- 'hs'
performance_df <- rbind(performance_df_hao, performance_df_mouse)
performance_df <- subset(performance_df, method_norm_combi %in% method_parameter_df$method_norm_combi)
performance_df <- performance_df %>%
  mutate(bulk_ds = recode(bulk_ds,
                          "chen" = "chen",
                          "chen-simulation" = "chenSim",
                          "chen-simulation-nobias" = "chenSim (no bias)",
                          "finotello" = "finotello",
                          "finotello-simulation" = "finotelloSim",
                          "finotello-simulation-nobias" = "finotelloSim (no bias)",
                          "hoek" = "hoek",
                          "hoek-simulation" = "hoekSim",
                          "hoek-simulation-nobias" = "hoekSim (no bias)",
                          "petitprez" = "petitprez",
                          "petitprez-simulation" = "petitprezSim",
                          "petitprez-simulation-nobias" = "petitprezSim (no bias)"))
performance_df <- performance_df %>% subset(bulk_ds %in% c('chenSim', 'chenSim (no bias)',
                                                           'finotelloSim', 'finotelloSim (no bias)',
                                                           'hoekSim', 'hoekSim (no bias)',
                                                           'petitprezSim', 'petitprezSim (no bias)'))
performance_df$withBias <- ifelse(grepl('no bias', performance_df$bulk_ds), F, T)

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

fractions_df_mouse <- get_fractions('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_mouse/', '0', 'tabula-muris', method_parameter_df)
fractions_df_mouse$org <- 'mm'
fractions_df_hao <- get_all_fractions('/nfs/data/omnideconv_benchmarking_clean/benchmark_results/results_main/', '0', 'hao-sampled-3', method_parameter_df)
fractions_df_hao$org <- 'hs'
fractions_df <- rbind(fractions_df_hao, fractions_df_mouse)
fractions_df <- subset(fractions_df, method_norm_combi %in% method_parameter_df$method_norm_combi)
fractions_df <- fractions_df %>% dplyr::rename(
  'cell_type' = 'celltype'
)
fractions_df <- fractions_df %>%
  mutate(bulk_ds = recode(bulk_ds,
                          "chen" = "chen",
                          "chen-simulation" = "chenSim",
                          "chen-simulation-nobias" = "chenSim (no bias)",
                          "finotello" = "finotello",
                          "finotello-simulation" = "finotelloSim",
                          "finotello-simulation-nobias" = "finotelloSim (no bias)",
                          "hoek" = "hoek",
                          "hoek-simulation" = "hoekSim",
                          "hoek-simulation-nobias" = "hoekSim (no bias)",
                          "petitprez" = "petitprez",
                          "petitprez-simulation" = "petitprezSim",
                          "petitprez-simulation-nobias" = "petitprezSim (no bias)"))
fractions_df <- fractions_df %>% subset(bulk_ds %in% c('chenSim', 'chenSim (no bias)',
                                                       'finotelloSim', 'finotelloSim (no bias)',
                                                       'hoekSim', 'hoekSim (no bias)',
                                                       'petitprezSim', 'petitprezSim (no bias)'))
fractions_df$withBias <- !grepl('no bias', fractions_df$bulk_ds)
fractions_df$base_dataset <- gsub(pattern = ' \\(no bias\\)', replacement = '', x = fractions_df$bulk_ds)
fractions_df$difference <- fractions_df$fraction.estimate - fractions_df$fraction.true
fractions_df$difference_norm <- (fractions_df$fraction.estimate - fractions_df$fraction.true) / fractions_df$fraction.true

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


#### plot figure S5a ####

performance_df$base_dataset <- gsub(pattern = ' \\(no bias\\)', replacement = '', x = performance_df$bulk_ds)
df <- performance_df %>% subset(cell_type == 'all')
df2 <- df[,.(RMSE_prop = RMSE[withBias] - RMSE[!withBias]), by =c('cell_type','method','base_dataset')]

fig_s5a <- ggplot(df2, aes(x=method, y=base_dataset, fill=RMSE_prop))+
  geom_tile()+
  coord_flip()+
  labs(fill='RMSE - RMSE_nobias')+
  scale_fill_gradient2(mid = "white", high = "#8a7495", low="#7f9574") +
  xlab('')+ylab('')+theme_minimal()

ggsave(filename = 'visualizations_final/supplement/fig_s5a.pdf', fig_s5a, height=4, width = 6)


#### plot figure S5b ####

result <- fractions_df %>%
  group_by(cell_type, method, base_dataset, org) %>%
  do({
    # Perform Wilcoxon test for each group
    group_data <- .
    
    # Split data based on 'withBias' levels
    delta_bias <- group_data$difference[group_data$withBias == TRUE]
    delta_nobias <- group_data$difference[group_data$withBias == FALSE]
    
    # Wilcoxon test between True and False levels of 'withBias'
    wilcox_result <- wilcox.test(
      x = abs(delta_nobias),
      y = abs(delta_bias),
      alternative = 'greater',
      paired = T
    )
    
    # Extract relevant information
    data.frame(
      cell_type = unique(group_data$cell_type),
      method = unique(group_data$method),
      W_statistic = wilcox_result$statistic,
      p_value = wilcox_result$p.value,
      logp = -log10(wilcox_result$p.value),
      avg_deltabias = mean(delta_bias),
      avg_deltanobias = mean(delta_nobias),
      avg_deltabias_abs = mean(abs(delta_bias)),
      avg_deltanobias_abs = mean(abs(delta_nobias)),
      med_deltabias = median(delta_bias),
      med_deltanobias = median(delta_nobias),
      med_deltabias_abs = median(abs(delta_bias)),
      med_deltanobias_abs = median(abs(delta_nobias))
    )
  })
result$padj <- p.adjust(result$p_value, method = 'fdr')
result$signif = ifelse(result$padj < 0.05, '*','')

fig_s5b <- ggplot(result, aes(x=cell_type, y=method, fill=med_deltabias))+
  geom_tile(color='black')+
  geom_text(aes(label = signif), size = 4)+
  facet_wrap(~base_dataset)+
  scale_fill_gradient2(mid = "white", high = "#8a7495", low="#7f9574") +
  #scale_fill_viridis_c()+
  labs(x = "", y = "", fill='Median delta_bias')+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_flip()

ggsave(filename = 'visualizations_final/supplement/fig_s5b.pdf', plot = fig_s5b, height=6, width = 7)
ggsave(filename = 'visualizations_final/supplement/fig_s5b.png', plot = fig_s5b, height=6, width = 7)
