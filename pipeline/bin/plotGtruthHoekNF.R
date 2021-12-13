#!/usr/bin/Rscript
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
results <- args[[2]]

library(ggplot2)
library(tidyverse)


remapCelltypesPG <- function(celltype){
  if(grepl(pattern = "T.cell", x = celltype)){
    return("T cell")
  } else if(grepl(pattern = "Tregs", x = celltype)){
    return("T cell")
  } else if(grepl(pattern = "Tcell", x = celltype)){
    return("T cell")
  } else if(grepl(pattern = "NK.cell", x = celltype)){
    return("NK")
  } else if(grepl(pattern = "B.cell", x = celltype)){
    return("B cell")
  } else if(grepl(pattern = "Bcells", x = celltype)){
    return("B cell")
  } else if(grepl(pattern = "Monocytes", x = celltype)){
    return("Monocyte")
  } else if(grepl(pattern = "mono", x = celltype)){
    return("Monocyte")
  } else if(grepl(pattern = "NK", x = celltype)){
    return("NK")
  } else if(grepl(pattern = "DC", x = celltype)){
    return("Dendritic cell")
  } else if(grepl(pattern = "Dendritic.cells", x = celltype)){
    return("Dendritic cell")
    #} else if(grepl(pattern = "Macro", x = celltype)){
    #  return("mono")
  } else {
    return("other")
  }
}


transformDataForPlottingPG <- function(deconv_df, reference, remap=FALSE, addAll=FALSE){
  deconv <- as.data.frame(deconv_df) %>% rownames_to_column(var="sample") %>% pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value")
  if(remap){
    deconv <- mutate(deconv, celltype= sapply(celltype, function(x) remapCelltypesPG(x))) %>% group_by(sample, celltype) %>% 
      summarise(predicted_value= sum(predicted_value))
  }
  ref <- reference %>% rownames_to_column(var="sample") %>% pivot_longer(!"sample", names_to = "celltype", values_to = "true_value")
  if(remap){
    ref <- mutate(ref, celltype= sapply(celltype, function(x) remapCelltypesPG(x)))
  }
  df <- inner_join(deconv, ref, by=c("celltype", "sample")) %>% mutate(facet = celltype)
  if(addAll){
     return(rbind(df, mutate(df, facet = "all")))
  } else {
    return(df)
  }
}

plotQuantiseqVsGtruthPG <- function(df, dataSetName){
  full <- df[! df$celltype %in% c("other", "Other"),] #%>% mutate(facet = celltype) %>% rbind( mutate(df, facet = "all"))
  facetLevels <- unique(full$facet)
  full$facet <- factor(full$facet, levels=facetLevels)#Bcells = "B cell", mono="Monocyte", NK="NK cell", Tcell="T cell", all="all")
  g <- ggplot(full, aes(x=true_value, y=predicted_value, color=celltype))+geom_point()+
    geom_abline()+ggpubr::stat_cor(label.sep = "\n", size=2.5, color="black", label.x.npc = 0.01)+
    xlab("Ground truth cell fractions")+ylab("quanTIseq cell fractions")
  ggsave(paste("/nfs/proj/omnideconv_benchmarking/quantiseqGtruth_", dataSetName, ".jpg", sep=""), g, width = 7.5, height = 5)
}

compareGroundTruth <- function(filepath, facs, addAll=FALSE){
	data <- NULL
  for(method in c("bisque", "cibersortx", "dwls", "momf", "scaden", "scdc")){
    decMatrix <- data.frame(get(paste("deconv_", method, sep="")))
    data <- transformDataForPlottingPG(decMatrix, facs, remap = TRUE, addAll = addAll) %>% mutate(method = method) %>% rbind(data)
  }
  p <- ggplot(data, aes(x=true_value, y=predicted_value, color=celltype))+geom_point()+geom_abline()+
    ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+facet_grid("facet~method")+
    theme(axis.text.x = element_text(angle=90), legend.position = "bottom")
  p
  ggsave(filename = paste(filepath, "comparisonGtruth.jpg", sep=""),p)
}

compareGroundTruth("/nfs/proj/omnideconv_benchmarking/hoek_all_", hoek_pbmc_facs, addAll = TRUE)
