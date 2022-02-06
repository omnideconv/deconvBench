#!/usr/bin/Rscript
"Usage: 
  plotSimulationCorrelationNF.R <results> <bulk> <remapping_sheet>
Options:
<results> list of deconvolution result rds files
<bulk> simulated bulk
<remapping_sheet> excel mapping sheet with remapping" -> doc

args <- docopt::docopt(doc)
deconvolution_results <- unlist(strsplit(gsub("\\[", "", gsub("\\]", "", args$results)), ", "))
bulk <- unlist(strsplit(gsub("\\[", "", gsub("\\]", "", args$bulk)), ", "))
remapping_sheet <- args$remapping_sheet
#load deconv results
# deconv_results <- list("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD4+_lambrechts_bisque.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD8+_lambrechts_bisque.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_B cell_lambrechts_bisque.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD4+_lambrechts_cibersortx.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD8+_lambrechts_cibersortx.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_B cell_lambrechts_cibersortx.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD4+_maynard_bisque.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD8+_maynard_bisque.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_B cell_maynard_bisque.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD4+_maynard_cibersortx.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_T cell CD8+_maynard_cibersortx.rds", 
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_B cell_maynard_cibersortx.rds")
# transformDataForPlotting <- function(deconv_df, reference, remap=FALSE, addAll=FALSE){
#   mapping <- remapCelltypesTree(facs_celltypes = rownames(reference), deconv_celltypes = colnames(deconv_df))
#   deconv <- as.data.frame(deconv_df) %>% tibble::rownames_to_column(var="sample") %>% pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value")
#   deconv <- merge(deconv, mapping, by.x = "celltype", by.y = "child_type")
#   deconv <- deconv %>% dplyr::group_by(sample, parent_type) %>%
#     dplyr::summarise(predicted_value= sum(predicted_value)) %>% dplyr::rename(celltype=parent_type)
#   ref <- reference %>% tibble::rownames_to_column(var="celltype") %>% pivot_longer(!"celltype", names_to = "sample", values_to = "true_value")
#   df <- dplyr::inner_join(deconv, ref, by=c("celltype", "sample")) %>% dplyr::mutate(facet = celltype)
#   if(addAll){
#     return(rbind(df, dplyr::mutate(df, facet = "all")))
#   } else {
#     return(df)
#   }
# }
# 
# 
library(tidyr)
library(dplyr)
library(ggplot2)

loadData <- function(deconv_results, bulk, tuples, allData=FALSE){
  data <- NULL
  for(res in deconv_results){
    result <- as.data.frame(readRDS(res))
    info <- unlist(strsplit(gsub(".rds", "", gsub("deconvolution_", "", lapply(res, basename))), "_"))
    facs <- as.data.frame(readRDS(tuples[tuples$scenario==info[1]&tuples$dataset==info[2], "path_bulk"][[1]])$cell_fractions) 
    if(info[2]=="maynard"){
      colnames(facs) <- c("Monocyte conventional", "Macrophage", "Endothelial cell", 
                          "Mast cell", "B cell", "T cell CD8+", "T cell CD4+", 
                          "Epithelial cell", "other cell", "NK cell",
                          "Monocyte non-conventional", "Myeloid dendritic cell", 
                          "Cancer associated fibroblast", "Club cell", "B cell plasma", 
                          "T cell regulatory (Tregs)", "Plasmacytoid dendritic cell", 
                          "Langerhans cell", "Goblet cell")
      
    } else {
      colnames(facs) <- c("B cell", "T cell regulatory (Tregs)", "other cell", 
                          "T cell CD8+", "Monocyte non-conventional", 
                          "Mast cell", "T cell CD4+", "Myeloid dendritic cell", 
                          "Cancer associated fibroblast", "Macrophage", "Club cell", 
                          "B cell plasma", "Goblet cell", "NK cell", "Langerhans cell", 
                          "Endothelial cell", "Monocyte conventional", 
                          "Plasmacytoid dendritic cell", "Epithelial cell")
    }
    facs <- facs %>% tibble::rownames_to_column(var="sample") %>%
      pivot_longer(!"sample", names_to="celltype", values_to = "true_value")
    result <- result %>% tibble::rownames_to_column(var="sample") %>%
      pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value") %>%
      mutate(dataset = info[2], method = info[3], combination = paste(info[2], info[3], sep="_"))
    full <- merge(result, facs, by=c("sample", "celltype")) %>% mutate(facet = celltype)
    data <- rbind(data, full)
    if(allData){
      data <- rbind(data, (full%>% mutate(facet="all")))
    }
  }
  return(data)
}


# ##remap celltype annotations of facs##
# source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
# celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
#                                                celltype_annotations = rownames(facs), 
#                                                method_ds = dataset_name)
# rownames(facs) <- celltype_annotations
# 
# #decMatrix <- data.frame(readRDS(resultVec[1]))
# #colnames(decMatrix) <- gsub("\\.", " ", colnames(decMatrix))
# ###code von fig bei 4.4.1 is gut, ausserdem heatmap in 3.4
# 
# 
# library(tidyr)
# library(ggplot2)
# 
# 
compareGroundTruth <- function(filepath, data, all=FALSE){ #TODO which facetting do we want?
  p <- ggplot(data, aes(x=true_value, y=predicted_value, color=celltype))+geom_point()+geom_abline()+
    ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+facet_grid("facet~combination")+
    theme(axis.text.x = element_text(angle=90))#, legend.position = "bottom")
  p
  if(all){
    ggsave(filename = file.path(filepath, "comparisonGtruth_all.jpeg"), p, height = 16, width = 16)
  } else {
    ggsave(filename = file.path(filepath, "comparisonGtruth.jpeg"), p, height = 16, width = 16)
  }
}

correlationPlots <- function(tuples, data, filepath){
  par(mfrow=c(length(unique(tuples$dataset)), length(unique(tuples$method))), mar=c(2,2,2,0))
  pdf(paste("singleCorrelationMatrices.pdf", sep=""))
  lapply(tuples$path_deconvolution, function(combination_path) {
    result <- as.data.frame(readRDS(combination_path))
    info <- unlist(strsplit(gsub(".rds", "", gsub("deconvolution_", "", lapply(combination_path, basename))), "_"))
    facs <- as.data.frame(readRDS(tuples[tuples$scenario==info[1]&tuples$dataset==info[2], "path_bulk"][[1]])$cell_fractions)
    if(info[2]=="maynard"){
      colnames(facs) <- c("Monocyte conventional", "Macrophage", "Endothelial cell", 
                          "Mast cell", "B cell", "T cell CD8+", "T cell CD4+", 
                          "Epithelial cell", "other cell", "NK cell",
                          "Monocyte non-conventional", "Myeloid dendritic cell", 
                          "Cancer associated fibroblast", "Club cell", "B cell plasma", 
                          "T cell regulatory (Tregs)", "Plasmacytoid dendritic cell", 
                          "Langerhans cell", "Goblet cell")
      
    } else {
      colnames(facs) <- c("B cell", "T cell regulatory (Tregs)", "other cell", 
                          "T cell CD8+", "Monocyte non-conventional", 
                          "Mast cell", "T cell CD4+", "Myeloid dendritic cell", 
                          "Cancer associated fibroblast", "Macrophage", "Club cell", 
                          "B cell plasma", "Goblet cell", "NK cell", "Langerhans cell", 
                          "Endothelial cell", "Monocyte conventional", 
                          "Plasmacytoid dendritic cell", "Epithelial cell")
    }
    
    corMatrix <- cor(result, facs)
    if(info[2]=="lambrechts"&info[3]=="bisque"){
      corrplot::corrplot(corMatrix,
                         is.corr = TRUE, tl.col = "black", #tl.pos = 'n', #order = 'hclust',
                         title = paste(info[2], " - ", info[3]), mar=c(0,0,2,0))
    } else {
      corrplot::corrplot(corMatrix,
                       is.corr = TRUE, tl.col = "black", #tl.pos = 'n', #order = 'alphabet',
                       title = paste(info[2], " - ", info[3]), mar=c(0,0,2,0))
    }
  })
  dev.off()
  
  subset <- data %>% select(-sample) %>% group_by(celltype, combination) %>% 
    mutate(cor = cor(true_value, predicted_value)) %>% 
    distinct(celltype, cor, combination)
  #matr <- subset %>% pivot_wider(names_from = combination, values_from = cor) %>% tibble::column_to_rownames("celltype")
  g <- ggplot(subset, aes(combination, celltype))+
    geom_tile(aes(fill=cor))+
    geom_text(aes(label=round(cor, 2)))+
    scale_fill_gradient(low = "white", high = "#1b98e0")+
    theme(axis.text.x = element_text(angle=90))
  ggsave(file.path(filepath, "celltypesCorrelationMatrix.jpeg"), g)
}

resTable <- tibble(path_deconvolution = deconvolution_results) %>%
  mutate(results = gsub(".rds", "", gsub("deconvolution_", "", lapply(path_deconvolution, basename)))) %>%
  separate(results, c("scenario", "dataset", "method"), sep="_")
bulkTable <- tibble(path_bulk = bulk) %>%
  mutate(results = gsub(".rds", "", gsub("simulatedBulk_", "", lapply(path_bulk, basename)))) %>%
  separate(results, c("scenario", "dataset"), sep="_") 
tuples <- merge(resTable, bulkTable, by=c("scenario", "dataset"))

data_all <- loadData(deconvolution_results, bulk, tuples, TRUE)
data <- loadData(deconvolution_results, bulk, tuples)
compareGroundTruth(filepath = ".", data = data_all, TRUE)
correlationPlots(tuples, data, ".")

# deconvolution_results <- c("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/work/c5/e4b043d91d81a55d941480746964dc/deconvolution_uniform_maynard_bisque.rds", 
#                            "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/work/d0/481cc89a02184a7cf8c5b7233a08d9/deconvolution_uniform_lambrechts_bisque.rds",
#                            "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/work/4c/e297275676c1c4d28af4d497f947b5/deconvolution_uniform_maynard_cibersortx.rds",
#                            "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/work/fb/bffe52785b09346ef689f7ef4d999f/deconvolution_uniform_lambrechts_cibersortx.rds")