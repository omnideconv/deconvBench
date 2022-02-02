#!/usr/bin/Rscript
"Usage: 
  plotSimulationCorrelationNF.R <results> <simulation>
Options:
<results> list of deconvolution result rds files
<facs> folder of matrix with facs fractions
<dataset_name> name of rnaseq dataset
<remapping_sheet> excel mapping sheet with remapping" -> doc

print(doc)
args <- docopt::docopt(doc)
results <- args$results
resultVec <- strsplit(gsub("\\[", "", gsub("\\]", "", results)), ", ")[[1]]
print(resultVec)
remapping_sheet <- args$remapping_sheet
dataset_name <- args$dataset_name
#TODO
#resave pbmc data as rds not RData!!
load(file.path(args$facs, dataset_name, paste(dataset_name, "_pbmc_facs.RData", sep="")))
facs <- get(paste(dataset_name, "_pbmc_facs", sep=""))
#rownames(facs) <- gsub("\\.", " ", rownames(facs))

##remap celltype annotations of facs##
source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                               celltype_annotations = rownames(facs), 
                                               method_ds = dataset_name)
rownames(facs) <- celltype_annotations

# library("ggplot2")
# g <- ggplot(data.frame(x=c(1,2,3), y=c(1,4,7)), aes(x, y))+geom_point()
# g
# ggsave("out.jpeg")

#decMatrix <- data.frame(readRDS(resultVec[1]))
#colnames(decMatrix) <- gsub("\\.", " ", colnames(decMatrix))
###code von fig bei 4.4.1 is gut, ausserdem heatmap in 3.4


library(tidyr)
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

compareGroundTruth <- function(filepath, facs, addAll=FALSE){
  data <- NULL
  for(i in 1:length(resultVec)){
    decMatrix <- data.frame(readRDS(resultVec[i]))
    colnames(decMatrix) <- gsub("\\.", " ", colnames(decMatrix))
    #bad fix: methods remove special characters so mapping gets screwed
    #this only works for finotello dataset now...
    # colnames(decMatrix) <- c("B cell", "B cell plasma", "Cancer associated fibroblast", 
    #                          "Club cell", "Endothelial cell", "Epithelial cell", "Goblet cell", 
    #                          "Langerhans cell", "Macrophage", "Mast cell", "Monocyte conventional", 
    #                          "Monocyte non-conventional", "Myeloid dendritic cell", 
    #                          "NK cell", "other cell", "Plasmacytoid dendritic cell", 
    #                          "T cell CD4+", "T cell CD8+", "T cell regulatory (Tregs)")
    # 
    combination <- gsub(".rds", "", gsub("deconvolution_", "", basename(resultVec[i])))
    print(combination)
    data <- transformDataForPlotting(deconv_df=decMatrix, reference=facs, remap = TRUE, addAll = addAll) %>% dplyr::mutate(combination = combination) %>% rbind(data)
  }
  p <- ggplot(data, aes(x=true_value, y=predicted_value, color=celltype))+geom_point()+geom_abline()+
    ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+facet_grid("facet~combination")+
    theme(axis.text.x = element_text(angle=90), legend.position = "bottom")
  p
  ggsave(filename = file.path(filepath, "comparisonGtruth.jpeg"),p)
}

compareGroundTruth(filepath = ".", facs = facs, addAll = TRUE)


### adjust this (see finalReport)
corMatrix <- cor(deconvolution, simulated_bulk$cell_fractions)
pdf(paste("correlationMatrix_", scenario, "_", method, "_", sc_ds, "_", 
          rnaseq_ds, ".pdf", sep=""))
corrplot::corrplot(corMatrix, 
                   is.corr = TRUE, tl.col = "black", order = 'alphabet', 
                   title = paste(scenario, "-", method))
dev.off()
