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

deconv_results <- list.files("/nfs/proj/omnideconv_benchmarking/pipelines/results_spillover", pattern = "deconvolution_spillover", full.names = TRUE)

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

getPalette = colorRampPalette(RColorBrewer::brewer.pal(11, "Paired"))
colors = getPalette(length(unique(data$true_celltype)))
names(colors) = unique(data$true_celltype)

signalRatiosBox <- function(resultDfIn, filename){
  signals = resultDfIn %>%
    mutate(signal_noise=if_else(celltype == true_celltype, "signal", "noise")) %>%
    group_by(dataset, rnaseq_type, sample, method, true_celltype, signal_noise) %>%
    summarise(estimate=sum(predicted_value)) %>%
    spread(signal_noise, estimate) %>%
    mutate(noise_ratio = noise/(noise+signal)) %>%
    mutate(signal_ratio = signal/(noise+signal)) %>%
    ungroup() %>%
    na.omit()
  
  ggplot(signals, aes(x=method, y=signal_ratio, fill=method)) +
    geom_boxplot(position = position_dodge()) +
    facet_grid(dataset+rnaseq_type~true_celltype, labeller = label_wrap_gen(width=10)) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "bottom") +
    labs(y="signal ratio")
  ggsave(filename, width = 19, height = 9)
}

### code adopted from immunedeconv benchmarking / spillover analysis
makeChordDiagrams <- function(resultDfIn, overviewTable, pdfName){
  par(mar=rep(0.5, 4))
  circos.par(cell.padding = rep(0, 4))
  pdf(pdfName, width = 20, height = 9)
  layout(matrix(seq(1, length(unique(overviewTable$dataset)) * length(unique(overviewTable$method)) * length(unique(overviewTable$rnaseq_type))), 
                length(unique(overviewTable$dataset)) * length(unique(overviewTable$rnaseq_type)),
                length(unique(overviewTable$method)) 
                ))
  lapply(unique(overviewTable$method), function(method){
    lapply(unique(overviewTable$dataset), function(dataset) {
      lapply(unique(overviewTable$rnaseq_type), function(type) {
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
          #group_by(method, celltype, true_celltype) %>%
          #summarise(estimate = mean(estimate)) %>%
          #group_by(method) %>%
          mutate(type = ifelse(celltype == true_celltype, "signal", "noise")) %>%
          group_by(method, type) %>%
          summarise(estimate = sum(estimate)) %>%
          spread(type, estimate) %>%
          mutate(noise_ratio = noise/(signal+noise), signal_ratio = signal/(signal+noise)) %>%
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
  # lgd_links = ComplexHeatmap::Legend(at = c(-2, -1, 0, 1, 2), 
  #                    title_position = "topcenter", title = "cell type")
  # ComplexHeatmap::draw(lgd_links, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = "centre")
  dev.off()
}

signalRatiosBox(resultDfIn = data, filename = "spillover_signalRatio.jpeg")
makeChordDiagrams(resultDfIn = data, overviewTable = resTable, pdfName = "spilloverChordDiagram.pdf")

x = data %>%
  #filter(method == "dwls" & rnaseq_type == "counts") %>%
  group_by(method, true_celltype, celltype) %>%
  summarise(frac = sum(predicted_value)/n())
ggplot(x, aes(x=true_celltype, y=frac)) +
  geom_bar(aes(fill=celltype), stat="identity") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~method)+labs(x="true celltype", y="prediction", fill="predicted celltype")
ggsave("spillover_fractionsPerMethodBar.jpeg", width=13, height = 8)
