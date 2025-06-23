library(tidyverse)
library(dplyr)
library(circlize)
library(ggpubr)
library(cowplot)

base_path <- '/vol/omnideconv_input/results_allen_resolution'
files <- list.files(base_path, full.names = TRUE, recursive = TRUE, pattern = 'deconvolution.rds')

aggr_map <- read.table('/vol/omnideconv_results/omnideconv_data/allen/cell_type_mappings.csv', sep=',', header=TRUE)
names(aggr_map) <- c('fine', 'normal')

extract_metadata <- function(path) {
  method <- str_extract(path, "(?<=results_allen_resolution/)[^/]+")
  replicate <- str_extract(path, "replicate_\\d+")
  annotation_type <- str_extract(path, "(?<=allen_resolution_analysis_sim_)[^/]+")
  tibble(path = path, method = method, replicate = replicate, annotation_type = annotation_type)
}

results_long <- map_dfr(files, function(file_path) {
  meta <- extract_metadata(file_path)
  obj <- readRDS(file_path)

  pred_mat <- obj[[1]]
  if (meta$method == "scaden" && meta$annotation_type == "fine_annot") {
    new_names <- colnames(pred_mat)
    new_names <- str_replace_all(new_names, c("L2.3.IT" = "L2/3 IT",
                                              "L5.6.NP" = "L5/6 NP",
                                              "Micro.PVM" = "Micro-PVM"))

    new_names <- str_replace_all(new_names, "\\.", " ")
    colnames(pred_mat) <- new_names
  }

  if (meta$method == "dwls" && meta$annotation_type == "fine_annot") {
    new_names <- colnames(pred_mat)
    new_names <- str_replace_all(new_names, c("L2-3_IT" = "L2/3 IT",
                                              "L5-6_NP" = "L5/6 NP"))

    new_names <- str_replace_all(new_names, "\\_", " ")
    colnames(pred_mat) <- new_names
  }

  predicted <- pred_mat %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    pivot_longer(-sample, names_to = "cell_type", values_to = "predicted_fraction")

  true_mat <- obj[[2]]
  if (meta$annotation_type == "normal_annot") {
    rownames(true_mat)[rownames(true_mat) == "Micro-PVM"] <- "Micro"
    true_mat <- rbind(true_mat, Endo = 1 - colSums(true_mat))
  }

  true <- true_mat %>%
    as.data.frame() %>%
    rownames_to_column("cell_type") %>%
    pivot_longer(-cell_type, names_to = "sample", values_to = "true_fraction")

  base <- predicted %>%
    inner_join(true, by = c("sample", "cell_type")) %>%
    mutate(method = meta$method,
           replicate = meta$replicate,
           annotation_type = meta$annotation_type,
           path = file_path)

  # If fine_annot, add aggregated fine_aggr
  if (meta$annotation_type == "fine_annot") {
    aggr <- base %>%
      inner_join(aggr_map, by = c("cell_type" = "fine")) %>%
      group_by(sample, normal, method, replicate, path) %>%
      summarise(predicted_fraction = sum(predicted_fraction),
                true_fraction = sum(true_fraction),
                .groups = "drop") %>%
      rename(cell_type = normal) %>%
      mutate(annotation_type = "fine_aggr")

    base <- bind_rows(base, aggr)
  }

  return(base)
})

glimpse(results_long)




library(dplyr)


results_long <- results_long %>%
  mutate(coarse_cell_type = case_when(annotation_type %in% c("normal_annot", "fine_aggr") ~ cell_type,
                                      annotation_type == "fine_annot" ~ aggr_map$normal[match(cell_type, aggr_map$fine)],
                                      TRUE ~ NA_character_))

metrics_summary <- results_long %>%
  group_by(annotation_type, method, cell_type, coarse_cell_type) %>%
  summarise(pearson_r = cor(predicted_fraction, true_fraction, method = "pearson"),
            rmse = sqrt(mean((predicted_fraction - true_fraction)^2)),
            .groups = "drop")

data_to_plot <- metrics_summary[metrics_summary$coarse_cell_type %in% c('Excit', 'Inhib'), ]

data_to_plot$annotation_type[data_to_plot$annotation_type == "fine_aggr"] <- "fine (aggregated)"
data_to_plot$annotation_type[data_to_plot$annotation_type == "fine_annot"] <- "fine"
data_to_plot$annotation_type[data_to_plot$annotation_type == "normal_annot"] <- "coarse"



allen_palette <- c('Astro' = '#88CCEE',
                   'Endo' = '#c2c0b3',
                   'Micro' = '#fffb6a',
                   'Inhib' = '#44AA99',
                   'Excit' = '#CC6677',
                   'Oligo' = '#882255',
                   'OPC' = '#AA4499')


method_palette <- c('AutogeneS' = '#88CCEE',
                    'Bisque' = '#117733',
                    'BayesPrism' = '#c1772c',
                    'CIBERSORTx' = '#44AA99',
                    'DWLS' = '#CC6677',
                    'MuSiC' = '#332288',
                    'Scaden' = '#AA4499',
                    'SCDC' = '#c12c2c')

data_to_plot$method <- recode(data_to_plot$method,
                              'autogenes' = 'AutogeneS',
                              'bayesprism' = 'BayesPrism',
                              'bisque' = 'Bisque',
                              'cibersortx' = 'CIBERSORTx',
                              'dwls' = 'DWLS',
                              'music' = 'MuSiC',
                              'scaden' = 'Scaden',
                              'scdc' = 'SCDC')

plot1 <- data_to_plot[data_to_plot$annotation_type != 'fine', ] %>%
  ggplot(., aes(x=rmse, y=pearson_r, color=method, shape=annotation_type, group=method)) +
  geom_point(size=3) +
  geom_line(linetype='dashed') +
  theme_bw() +
  facet_wrap(.~cell_type) +
  scale_color_manual(values=method_palette) +
  labs(x='RMSE', y='R', shape="annotation type") +
  geom_hline(yintercept = 0.8, linetype = 'dashed') +
  geom_vline(xintercept = 0.1, linetype = 'dashed') +
  theme(legend.position = 'bottom')


plot2 <- data_to_plot[data_to_plot$annotation_type == 'fine', ] %>%
  ggplot(., aes(x=rmse, y=pearson_r, color=coarse_cell_type)) +
  geom_point(size=2) +
  theme_bw() +
  facet_wrap(.~method, ncol=4) +
  scale_color_manual(values=allen_palette) +
  labs(x='RMSE', y='R', color="coarse cell type") +
  geom_hline(yintercept = 0.8, linetype = 'dashed') +
  geom_vline(xintercept = 0.05, linetype = 'dashed') +
  rotate_x_text(angle=60) +
  theme(legend.position = 'bottom')




library(patchwork)


combined_plot <- plot2 / plot1 +
  plot_layout(heights = c(0.6, 0.4)) +
  plot_annotation(tag_levels = 'A')

ggsave(combined_plot, filename = './visualisation/supplement/Fig_s6.pdf', width = 13, height = 8)
