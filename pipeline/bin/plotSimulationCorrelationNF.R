#!/usr/bin/Rscript
"Usage: 
  plotSimulationCorrelationNF.R <results> <bulk> <remapping_sheet> [<type>]
Options:
<results> list of deconvolution result rds files
<bulk> simulated bulk
<remapping_sheet> excel mapping sheet with remapping
<type> simulation scenario" -> doc

args <- docopt::docopt(doc)
deconvolution_results <- unlist(strsplit(gsub("\\[", "", gsub("\\]", "", args$results)), ", "))
bulk <- unlist(strsplit(gsub("\\[", "", gsub("\\]", "", args$bulk)), ", "))
remapping_sheet <- args$remapping_sheet

#for testing
deconvolution_results <- list.files("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results", pattern = "deconvolution_truefraction", full.names = TRUE)
bulk <- list.files("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results", pattern = "simulatedBulk_true", full.names=TRUE)
r_threshold = 0.8
lower_label <- paste("R lower than/equal to", r_threshold)
higher_label <- "high"

library(tidyr)
library(dplyr)
library(ggplot2)
source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/plottingHelperNF.R")


loadData <- function(deconv_results, bulk, tuples, allData=FALSE){
  data <- NULL
  for(res in deconv_results){
    result <- as.data.frame(readRDS(res))
    if(length(result)>0){ ##results empty for eg scaden counts
      info <- unlist(strsplit(gsub(".rds", "", gsub("deconvolution_", "", basename(res))), "_"))
      facs <- as.data.frame(readRDS(tuples[tuples$scenario==info[1]&tuples$scset==info[2], "path_bulk"][[1]])$cell_fractions) 
      facs <- facs %>% tibble::rownames_to_column(var="sample") %>%
        pivot_longer(!"sample", names_to="celltype", values_to = "true_value")
      result <- result %>% tibble::rownames_to_column(var="sample") %>%
        pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value") %>%
        mutate(scset = info[2], sctype = info[3], method = info[4], combination = paste(info[2], info[3], info[4], sep="_"), 
               combinationForMatrix = paste(info[2], info[3], sep="_"))
      full <- merge(result, facs, by=c("sample", "celltype")) %>% mutate(facet = celltype)
      data <- rbind(data, full)
      if(allData){
        data <- rbind(data, (full%>% mutate(facet="all")))
      }
    }
  }
  return(data)
}

correlationPlots <- function(tuples, data, filepath){
  #par(mfrow=c(length(unique(tuples$scset)), length(unique(tuples$method))), mar=c(2,2,2,0))
  pdf(paste("comparisonGtruth_truefractions_singleCorrelationMatrices.pdf", sep=""))
  lapply(tuples$path_deconvolution, function(combination_path) {
    result <- as.data.frame(readRDS(combination_path))
    if(length(result)!=0){
      info <- unlist(strsplit(gsub(".rds", "", gsub("deconvolution_", "", basename(combination_path))), "_"))
      facs <- as.data.frame(readRDS(tuples[tuples$scenario==info[1]&tuples$scset==info[2], "path_bulk"][[1]])$cell_fractions)
      # if(info[2]=="maynard"){
      #   colnames(facs) <- c("Monocyte conventional", "Macrophage", "Endothelial cell", 
      #                       "Mast cell", "B cell", "T cell CD8+", "T cell CD4+", 
      #                       "Epithelial cell", "other cell", "NK cell",
      #                       "Monocyte non-conventional", "Myeloid dendritic cell", 
      #                       "Cancer associated fibroblast", "Club cell", "B cell plasma", 
      #                       "T cell regulatory (Tregs)", "Plasmacytoid dendritic cell", 
      #                       "Langerhans cell", "Goblet cell")
      #   
      # } else {
      #   colnames(facs) <- c("B cell", "T cell regulatory (Tregs)", "other cell", 
      #                       "T cell CD8+", "Monocyte non-conventional", 
      #                       "Mast cell", "T cell CD4+", "Myeloid dendritic cell", 
      #                       "Cancer associated fibroblast", "Macrophage", "Club cell", 
      #                       "B cell plasma", "Goblet cell", "NK cell", "Langerhans cell", 
      #                       "Endothelial cell", "Monocyte conventional", 
      #                       "Plasmacytoid dendritic cell", "Epithelial cell")
      # }
      
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
    }
  })
  dev.off()
  
  subset <- data %>% select(-sample) %>% group_by(celltype, combination) %>% 
    mutate(cor = cor(true_value, predicted_value)) %>% 
    distinct(celltype, cor, combinationForMatrix, method, combination)
  #matr <- subset %>% pivot_wider(names_from = combination, values_from = cor) %>% tibble::column_to_rownames("celltype")
  g <- ggplot(subset, aes(combination, celltype))+
    geom_tile(aes(fill=cor))+
    geom_text(aes(label=round(cor, 2)))+
    scale_fill_gradient(low = "white", high = "#1b98e0")+
    theme(axis.text.x = element_text(angle=90)) + 
    labs(x="combination", y="cell type")
  g
  ggsave(file.path(filepath, "comparisonGtruth_truefractions_celltypesCorrelationMatrix.jpeg"), g, width = 18, height=13)
  g <- ggplot(subset, aes(combinationForMatrix, celltype))+
    geom_tile(aes(fill=cor))+
    geom_text(aes(label=round(cor, 2)))+
    scale_fill_gradient(low = "white", high = "#1b98e0")+
    theme(axis.text.x = element_text(angle=90)) + 
    labs(x="combination", y="cell type") + 
    facet_wrap(~method, ncol = length(unique(data$method)))
  g
  ggsave(file.path(filepath, "comparisonGtruth_truefractions_celltypesCorrelationMatrix_facetMethods.jpeg"), g, width = 18, height=13)
}

resTable <- tibble(path_deconvolution = deconvolution_results) %>%
  mutate(results = gsub(".rds", "", gsub("deconvolution_", "", basename(path_deconvolution)))) %>%
  separate(results, c("scenario", "scset", "sctype", "method"), sep="_")
bulkTable <- tibble(path_bulk = bulk) %>%
  mutate(results = gsub(".rds", "", gsub("simulatedBulk_", "", basename(path_bulk)))) %>%
  separate(results, c("scenario", "scset"), sep="_") 
tuples <- merge(resTable, bulkTable, by=c("scenario", "scset"))

data_all <- loadData(deconv_results = deconvolution_results, bulk = bulk, tuples = tuples, TRUE)
data <- loadData(deconv_results = deconvolution_results, bulk = bulk, tuples = tuples)

#between samples
separateByCondition(filepath = ".", data = data, 
                    subsetVectorFull = data$scset, 
                    basename = "comparisonGtruth_truefractions_betweenSampleComparison_colorCelltype", 
                    shape = "sctype", color = "celltype")
separateByCondition(filepath = ".", data = data, 
                    subsetVectorFull = data$scset, 
                    basename = "comparisonGtruth_truefractions_betweenSampleComparison_colorCelltype_lm_freeScales", 
                    shape = "sctype", color = "celltype", addLM = TRUE, scales = "free")
separateByCondition(filepath = ".", data = data_all, 
                    subsetVectorFull = data_all$scset, 
                    basename = "comparisonGtruth_truefractions_betweenSampleComparison_colorCelltype_lm_freeScales_plusAllCondition", 
                    shape = "sctype", color = "celltype", addLM = TRUE, scales = "free")

#between scsets
separateByCondition(filepath = ".", data = data, 
                    subsetVectorFull = data$sctype, 
                    basename = "comparisonGtruth_truefractions_betweenScsets_lm", 
                    shape = "scset", addLM = TRUE, scales = "free")

compareGroundTruthAllCombinations(filepath = ".", data = data, 
                                  basename = "comparisonGtruth_truefractions_allCombinations", 
                                  shape = "sctype",
                                  scenario = "simulation")

allTogether(filepath = ".", data = data, groupingVar1 = "method", 
            groupingVar2 = "sctype", basename = "comparisonGtruth_truefractions_all", 
            shape = "scset")

allTogether(filepath = ".", data = data_all, groupingVar1 = "method", 
            groupingVar2 = "facet", basename = "comparisonGtruth_truefractions_all", 
            shape = "scset", width = 18, height=8, addLM = TRUE)

allTogether(filepath = ".", data = mutate(data, noFacet="all"), groupingVar1 = "noFacet", 
            groupingVar2 = "method", basename = "comparisonGtruth_truefractions_all_noFacet", 
            shape = "scset")
allTogetherOneFacetNoShape(filepath = ".", data = data, groupingVar1 = "method", 
                           basename = "comparisonGtruth_truefractions_OnlyMethod")
correlationPlots(tuples, data, ".")
