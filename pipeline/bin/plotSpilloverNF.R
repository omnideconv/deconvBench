#!/usr/bin/Rscript

"Usage: 
  plotSpilloverNF.R <deconv_results>
Options:
<deconv_results> list with spillover simulation deconv results" -> doc

args <- docopt::docopt(doc)

library(tidyr)
library(dplyr)
library(circlize)
escapeCelltypes <- function(celltype){
  return(gsub(" ", "x.x", celltype))
}
reEscapeCelltypes <- function(celltype){
  return(gsub("x.x", " ", celltype))
}
deconv_results <- unlist(strsplit(gsub("\\[", "", gsub("\\]", "", args$deconv_results)), ", "))

deconv_results <- list.files("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results", pattern = "deconvolution_spillover", full.names = TRUE)
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
#                        "/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_B cell_maynard_cibersortx.rds")#, 
#"/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results/deconvolution_spillover_maynard_autogenes.rds")
resTable <- tibble(path = deconv_results) %>% 
  mutate(results = gsub(".rds", "", gsub("deconvolution_", "", basename(path)))) %>% 
  separate(results, c("scenario", "dataset", "celltype", "rnaseq_type", "method"), sep="_") %>% 
  mutate(celltype=reEscapeCelltypes(celltype))
data <- NULL
for(res in deconv_results){
  result <- as.data.frame(readRDS(res))
  info <- unlist(strsplit(gsub(".rds", "", gsub("deconvolution_", "", basename(res))), "_"))
  result <- result %>% tibble::rownames_to_column(var="sample") %>% 
    pivot_longer(!"sample", names_to="celltype", values_to = "predicted_value") %>% 
    mutate(dataset = info[2], rnaseq_type = info[4], method = info[5]) %>% 
    mutate(true_value = ifelse(celltype == reEscapeCelltypes(info[3]), 1, 0), true_celltype = reEscapeCelltypes(info[3])) #merge with gold standard
  data <- rbind(data, result)
} 

getPalette = colorRampPalette(brewer.pal(11, "Paired"))
colors = getPalette(length(unique(data$true_celltype)))
names(colors) = unique(data$true_celltype)
### code adopted from immunedeconv benchmarking / spillover analysis
makeChordDiagrams <- function(resultDfIn, overviewTable, pdfName){
  # layout(matrix(seq(1, length(unique(overviewTable$dataset)) * length(unique(overviewTable$method)) * length(unique(overviewTable$rnaseq_type))), 
  #               length(unique(overviewTable$method)), 
  #               length(unique(overviewTable$dataset)) * length(unique(overviewTable$rnaseq_type))))
  par(mar=rep(0.5, 4))
  circos.par(cell.padding = rep(0, 4))
  pdf(pdfName, width = 20, height = 20)
  layout(matrix(seq(1, length(unique(overviewTable$dataset)) * length(unique(overviewTable$method)) * length(unique(overviewTable$rnaseq_type))), 
                length(unique(overviewTable$method)), 
                length(unique(overviewTable$dataset)) * length(unique(overviewTable$rnaseq_type))))
  lapply(unique(overviewTable$rnaseq_type), function(type){
    lapply(unique(overviewTable$dataset), function(dataset) {
      lapply(unique(overviewTable$method), function(method) {
        resultDf <- resultDfIn[resultDfIn$rnaseq_type==type,]
        migration = resultDf %>%
          filter(method == !!method, dataset == !!dataset) %>%
          group_by(method, true_celltype, celltype) %>%
          summarise(estimate = mean(predicted_value)) %>%
          ungroup()
        migration_mat = migration %>%
          select(-method) %>%
          spread(celltype, estimate) %>%
          as.data.frame() %>%
          tibble::column_to_rownames("true_celltype") %>%
          as.matrix()
        noise_ratio = migration %>%
          group_by(method, celltype, true_celltype) %>%
          summarise(estimate = mean(estimate)) %>%
          group_by(method) %>%
          mutate(type = ifelse(celltype == true_celltype, "signal", "noise")) %>%
          group_by(method, type) %>%
          summarise(estimate = sum(estimate)) %>%
          spread(type, estimate) %>%
          mutate(noise_ratio = noise/(signal+noise)) %>%
          ungroup()
        chordDiagram(migration_mat, directional = TRUE, transparency = .5,
                     grid.col = colors
                     #annotationTrack = c("track", "grid")
        )
        text(0, 0, method, cex = 2.3)
        text(0, -0.2, as.character(round(filter(noise_ratio, method == !!method) %>% pull(noise_ratio), 2)), cex=2)
        # first method, add title.
        if(method == unique(overviewTable$method)[1]) {
          title(paste(dataset, type, sep=" - "))
        }
      })
    })
  })
  dev.off()
}

makeChordDiagrams(resultDfIn = data, overviewTable = resTable, pdfName = "spilloverChordDiagram.pdf")
