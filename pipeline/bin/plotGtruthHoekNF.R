#!/usr/bin/Rscript
"Usage: 
  plotGtruthHoekNF.R <results> <facs> <remapping_sheet>
Options:
<results> list of deconvolution result rds files
<facs> folder of matrix with facs fractions
<remapping_sheet> excel mapping sheet with remapping" -> doc

print(doc)
args <- docopt::docopt(doc)
results <- args$results
resultVec <- strsplit(gsub("\\[", "", gsub("\\]", "", results)), ", ")[[1]]
print(resultVec)
remapping_sheet <- args$remapping_sheet

##remap celltype annotations of facs##
source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")

combinations <- lapply(resultVec, function(x){gsub(".rds", "", gsub("deconvolution_", "", basename(x)))})
rnaseqDatasets <- unique(lapply(combinations, function(x){unlist(strsplit(x, split="_"))[3]}))
facsList <- rnaseqDatasets
names(facsList) <- rnaseqDatasets
facsList <- lapply(facsList, function(x){
  load(file.path(args$facs, x, paste(x, "_pbmc_facs.RData", sep="")))
  facs <- get(paste(x, "_pbmc_facs", sep=""))
  celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                 celltype_annotations = rownames(facs), 
                                                 method_ds = x)
  rownames(facs) <- celltype_annotations
  facs
})

library(tidyr)
library(dplyr)
library(ggplot2)

transformDataForPlotting <- function(deconv_df, reference, remap=FALSE, addAll=FALSE){
  mapping <- remapCelltypesTree(facs_celltypes = rownames(reference), deconv_celltypes = colnames(deconv_df))
  deconv <- as.data.frame(deconv_df) %>% tibble::rownames_to_column(var="sample") %>% pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value")
  deconv <- merge(deconv, mapping, by.x = "celltype", by.y = "child_type")
  deconv <- deconv %>% dplyr::group_by(sample, parent_type) %>%
      dplyr::summarise(predicted_value= sum(predicted_value)) %>% dplyr::rename(celltype=parent_type)
  ref <- reference %>% tibble::rownames_to_column(var="celltype") %>% pivot_longer(!"celltype", names_to = "sample", values_to = "true_value")
  df <- dplyr::inner_join(deconv, ref, by=c("celltype", "sample")) %>% dplyr::mutate(facet = celltype)
  if(addAll){
    return(rbind(df, dplyr::mutate(df, facet = "all")))
  } else {
    return(df)
  }
}

loadData <- function(resultVector, facsListInput, addAll=FALSE){
  data <- NULL
  for(i in 1:length(resultVector)){
    decMatrix <- as.data.frame(readRDS(resultVec[i]))
    combination <- gsub(".rds", "", gsub("deconvolution_", "", basename(resultVector[i])))
    info <- unlist(strsplit(combination, split = "_"))
    method <- info[1]
    scset <- info[2]
    rnaset <- info[3]
    facsDf <- facsListInput[[rnaset]]
    data <- transformDataForPlotting(deconv_df=decMatrix, reference=facsDf, remap = TRUE, addAll = addAll) %>% 
      dplyr::mutate(combination = combination, method=method, scset=scset, bulk=rnaset) %>% 
      rbind(data)
  }
  return(data)
}

##TODO function within sample (for each sample one cor value - all celltypes all sc each method)

betweenSample <- function(filepath, data){
  for(scset in unique(data$scset)){
    p <- ggplot(data[data$scset==scset,], aes(x=true_value, y=predicted_value, color=celltype, shape = bulk))+
      geom_point()+geom_abline()+
      ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+
      facet_grid("facet~method")+
      theme(axis.text.x = element_text(angle=90), legend.position = "bottom")+
      labs(title = scset)
    p
    ggsave(filename = file.path(filepath, paste("comparisonGtruth_regularDeconv_betweenSampleComparison_", scset, ".jpeg", sep="")),p)
  }
}

compareGroundTruthOneRna <- function(filepath, data, addAll=FALSE){
  for(rnaSeqName in unique(data$bulk)){
    p <- ggplot(data[data$bulk == rnaSeqName,], aes(x=true_value, y=predicted_value))+
      geom_point(aes(color=celltype, shape=scset))+geom_abline()+
      ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+
      facet_grid("facet~method")+
      theme(axis.text.x = element_text(angle=90), legend.position = "bottom")+
      labs(title = rnaSeqName)
    p
    ggsave(filename = file.path(filepath, 
                                paste("comparisonGtruth_regularDeconv_betweenScsets_", rnaSeqName, ".jpeg", sep="")),p)
  }
  p <- ggplot(data, aes(x=true_value, y=predicted_value))+
    geom_point(aes(color=celltype, shape=bulk))+geom_abline()+
    ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+
    facet_grid("scset~method")+
    theme(axis.text.x = element_text(angle=90), legend.position = "bottom")
  p
  ggsave(filename = file.path(filepath, "comparisonGtruth_regularDeconv_betweenScsets_scatter.jpeg"),p)
  
  subset <- data %>% select(-sample) %>% group_by(scset, method) %>% 
    mutate(cor = cor(true_value, predicted_value)) %>% 
    distinct(celltype, cor, scset, method)
  #boxplot of pearson distribution
  g <- ggplot(subset, aes(cor, method))+
    geom_boxplot()+
    labs(x="pearson")+
    facet_wrap(~ celltype)
  g
  ggsave(filename = file.path(filepath, "comparisonGtruth_regularDeconv_betweenScsets_box.jpeg"),g)
  
  g <- ggplot(subset, aes(method, scset))+
    geom_tile(aes(fill=cor))+
    geom_text(aes(label=round(cor, 2)))+
    scale_fill_gradient(low = "white", high = "#1b98e0")+
    theme(axis.text.x = element_text(angle=90))
  g
  ggsave(file.path(filepath, "comparisonGtruth_regularDeconv_betweenScsets_correlation.jpeg"), g)
}

compareGroundTruthAllRna <- function(filepath, data, addAll=FALSE){
  p <- ggplot(data, aes(x=true_value, y=predicted_value, color=celltype))+
    geom_point()+geom_abline()+
    ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+
    facet_grid("facet~combination")+
    theme(axis.text.x = element_text(angle=90), legend.position = "bottom")
  p
  ggsave(filename = file.path(filepath, "comparisonGtruth_regularDeconv_all.jpeg"),p)
}

inputData <- loadData(resultVector = resultVec, facsListInput = facsList)
betweenSample(".", data = inputData)
compareGroundTruthOneRna(".", data=inputData)
compareGroundTruthAllRna(filepath = ".", data=inputData, addAll = TRUE)
