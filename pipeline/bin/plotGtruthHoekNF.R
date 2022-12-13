#!/usr/bin/Rscript
"Usage: 
  plotGtruthHoekNF.R <results> <facs> <remapping_sheet>
Options:
<results> list of deconvolution result rds files
<bulk_dataset> bulk dataset deconvolved
<norm_bulk_dataset> normalization of the bulk dataset
<facs> folder of matrix with facs fractions
<remapping_sheet> excel mapping sheet with remapping" -> doc

print(doc)
args <- docopt::docopt(doc)
results <- args$results
resultVec <- strsplit(gsub("\\[", "", gsub("\\]", "", results)), ", ")[[1]]
print(resultVec)
remapping_sheet <- args$remapping_sheet

#for testing
resultVec <- list.files("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/", 
                        pattern = "^deconvolution_autogenes|^deconvolution_bisque|^deconvolution_cdseq|^deconvolution_cibersortx|^deconvolution_cpm|^deconvolution_dwls|^deconvolution_momf|^deconvolution_music_|^deconvolution_scaden_", 
                        full.names=TRUE)

##remap celltype annotations of facs##
source("/vol/spool/bin/remapCelltypesNF.R")
source("/vol/spool/bin/plottingHelperNF.R")

combinations <- lapply(resultVec, function(x){gsub(".rds", "", gsub("deconvolution_", "", basename(x)))})
rnaseqDatasets <- unique(lapply(combinations, function(x){unlist(strsplit(x, split="_"))[4]}))
facsList <- rnaseqDatasets
names(facsList) <- rnaseqDatasets
facsList <- lapply(facsList, function(x){
  facs <- readRDS(file.path(args$facs, x, paste(x, "_facs.rds", sep="")))
  celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                 celltype_annotations = rownames(facs), 
                                                 method_ds = x)
  rownames(facs) <- celltype_annotations
  facs
})

#set threshold of R value!
r_threshold = 0.6
lower_label <- paste("R lower than/equal to", r_threshold)
higher_label <- "high"

library(tidyr)
library(dplyr)
library(ggplot2)

transformDataForPlotting <- function(deconv_df, reference, remap=FALSE, addAll=FALSE){
  mapping <- remapCelltypesTree(facs_celltypes = rownames(reference), deconv_celltypes = colnames(deconv_df))
  deconv <- as.data.frame(deconv_df) %>% tibble::rownames_to_column(var="sample") %>% 
    pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value")
  deconv <- merge(deconv, mapping, by.x = "celltype", by.y = "child_type")
  deconv <- deconv %>% dplyr::group_by(sample, parent_type) %>%
      dplyr::summarise(predicted_value= sum(predicted_value)) %>% 
    dplyr::rename(celltype=parent_type)
  ref <- as.data.frame(reference) %>% tibble::rownames_to_column(var="celltype") %>% 
    pivot_longer(!"celltype", names_to = "sample", values_to = "true_value")
  df <- dplyr::inner_join(deconv, ref, by=c("celltype", "sample")) %>% 
    dplyr::mutate(facet = celltype)
  if(addAll){
    return(rbind(df, dplyr::mutate(df, facet = "all")))
  } else {
    return(df)
  }
}

loadData <- function(resultVector, facsListInput, addAll=FALSE){
  data <- NULL
  for(i in 1:length(resultVector)){
    print(resultVector[i])
    decMatrix <- as.data.frame(readRDS(resultVector[i]))
    combination <- gsub(".rds", "", gsub("deconvolution_", "", basename(resultVector[i])))
    info <- unlist(strsplit(combination, split = "_"))
    method <- info[1]
    scset <- info[2]
    sctype <- info[3]
    rnaset <- info[4]
    facsDf <- facsListInput[[rnaset]]
    if(length(decMatrix)==0){
      #data <- rbind(data, data.frame(sample=NA, celltype=NA, predicted_value=NA, true_value=NA, facet=NA, combination = combination, method=method, scset=scset, 
      #                               sctype=sctype, bulk=rnaset))
    } else {
      data <- transformDataForPlotting(deconv_df=decMatrix, reference=facsDf, 
                                       remap = TRUE, addAll = addAll) %>% 
        dplyr::mutate(combination = combination, method=method, scset=scset, 
                      sctype=sctype, bulk=rnaset) %>% 
        rbind(data)
    }
  }
  return(data)
}


inputData <- loadData(resultVector = resultVec, facsListInput = facsList)
inputDataAll <- loadData(resultVector = resultVec, facsListInput = facsList, addAll = TRUE)

#between samples
separateByCondition(filepath = ".", data = inputDataAll, 
                    subsetVectorFull = inputDataAll$scset, 
                    basename = "comparisonGtruth_regularDeconv_betweenSampleComparison_plusAllCelltypes_colorType", 
                    shape = "bulk", color = "sctype")
separateByCondition(filepath = ".", data = inputData, 
                    subsetVectorFull = inputData$scset, 
                    basename = "comparisonGtruth_regularDeconv_betweenSampleComparison_colorType", 
                    shape = "bulk", color = "sctype")
separateByCondition(filepath = ".", data = inputData, 
                    subsetVectorFull = inputData$scset, 
                    basename = "comparisonGtruth_regularDeconv_betweenSampleComparison_colorType_lm_freeScales", 
                    shape = "bulk", color = "sctype", addLM = TRUE, scales = "free")
separateByCondition(filepath = ".", data = inputData, 
                    subsetVectorFull = inputData$scset, 
                    basename = "comparisonGtruth_regularDeconv_betweenSampleComparison_colorCelltype", 
                    shape = "bulk")
separateByCondition(filepath = ".", data = inputDataAll, 
                    subsetVectorFull = inputDataAll$scset, 
                    basename = "comparisonGtruth_regularDeconv_betweenSampleComparison_plusAllCelltypes_colorCelltype", 
                    shape = "bulk")

#between samples 

#between scsets
separateByCondition(filepath = ".", data = inputData, 
                    subsetVectorFull = inputData$bulk, 
                    basename = "comparisonGtruth_regularDeconv_betweenScsets_lm", 
                    shape = "scset", color = "sctype", scales = "free")

#healthy vs tumor
allTogether(filepath = ".", 
            data = subset(mutate(inputData, tumor = ifelse(grepl("vanderbilt", bulk), "tumor", "healthy")), grepl("T cell C", celltype)), 
            groupingVar1 = "celltype", groupingVar2 = "tumor",
            basename = "comparisonGtruth_regularDeconv_betweenHealthyTumor", 
            shape = "bulk", color = "bulk")
allTogether(filepath = ".", 
            data = subset(mutate(inputData, tumor = ifelse(grepl("vanderbilt", bulk), "tumor", "healthy")), grepl("T cell C", celltype)), 
            groupingVar1 = "method", groupingVar2 = "tumor",
            basename = "comparisonGtruth_regularDeconv_betweenHealthyTumor", 
            shape = "bulk", color = "celltype", addLM = TRUE)
allTogetherOneFacetNoShape(filepath = ".", 
                           data = subset(mutate(inputData, tumor = ifelse(grepl("vanderbilt", bulk), "tumor", "healthy")), grepl("T cell C", celltype)), 
                           groupingVar1 = "tumor", 
                           basename = "comparisonGtruth_regularDeconv_betweenHealthyTumor")


compareGroundTruthAllCombinations(filepath = ".", data = inputData, 
                                  basename = "comparisonGtruth_regularDeconv", shape = "sctype")

allTogether(filepath = ".", data = inputData, groupingVar1 = "celltype", 
            groupingVar2 = "combination", 
            basename = "comparisonGtruth_regularDeconv_all", shape = "scset")

allTogetherOneFacetNoShape(filepath = ".", data = inputData, 
                           groupingVar1 = "method", 
                           basename = "comparisonGtruth_regularDeconv_OnlyMethod")


