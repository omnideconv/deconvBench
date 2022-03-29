#!/usr/bin/Rscript
"Usage: 
  plotMinimumFractionNF.R <results> <remapping_sheet>
Options:
<results> list of deconvolution result rds files
<remapping_sheet> excel mapping sheet with remapping" -> doc

args <- docopt::docopt(doc)
deconvolution_results <- unlist(strsplit(gsub("\\[", "", gsub("\\]", "", args$results)), ", "))
#bulk <- unlist(strsplit(gsub("\\[", "", gsub("\\]", "", args$bulk)), ", "))
remapping_sheet <- args$remapping_sheet

escapeCelltypes <- function(celltype){
  return(gsub(" ", "x.x", celltype))
}
reEscapeCelltypes <- function(celltype){
  return(gsub("x.x", " ", celltype))
}

#for testing
deconvolution_results <- list.files("/nfs/proj/omnideconv_benchmarking/pipelines/results", pattern = "deconvolution_spikein.*Bx.xcell_.*tpm", full.names = TRUE)
#bulk <- list.files("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/results", pattern = "simulatedBulk_true", full.names=TRUE)
r_threshold = 0.8
lower_label <- paste("R lower than/equal to", r_threshold)
higher_label <- "high"

library(tidyr)
library(dplyr)
library(ggplot2)
source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/plottingHelperNF.R")


loadData <- function(deconv_results){
  data <- NULL
  bg_list <- list()
  mse_mean <- NULL
  celltypeOld <- "new"
  percentOld = 0
  sums <- NULL
  #for(res in deconv_results[-1]){ #1-5, 21-25, 31-35
  #for(i in c(1,2,3,4,5,21,22,23,24,25,31,32,33,34,35)){
  for(i in 1:15){
    res <- deconv_results[-1][i]
    print(res)
    result <- as.data.frame(readRDS(res))
    info <- unlist(strsplit(gsub(".rds", "", gsub("deconvolution_", "", basename(res))), "_"))
    celltype <- reEscapeCelltypes(info[4])
    # if(info[3]=="backgroundFor"){
    #   real <- data.frame(background = result[celltype], celltype=celltype, 
    #                      type=info[5], method=info[6])
    #   colnames(real) <- c("background", "celltype", "type", "method")
    #   bg_list[[paste(info[4], info[5], info[6], sep ="_")]] <- real
    # } else { #has 2 fields more
      percent <- as.integer(gsub("nrCells", "", info[5]))/10
      tmp <- data.frame(background = result[celltype], celltype=celltype, 
                        percent = percent, type=info[7], method=info[8])
      colnames(tmp) <- c("fraction", "celltype", "percent")
      print(summary(tmp$fraction))
      #bg_real <- bg_list[[paste(info[4], info[7], info[8], sep ="_")]]
      tmp$MSE <- (bg_real$background-tmp$fraction)^2
      #data <- rbind(data, tmp) 
      if(celltypeOld=="new"){
        celltypeOld = celltype
        print("changed celltype and percent")
        percentOld = percent
        sumsVec <- c()
        data <- tmp
        percentOld = percent
      } else if(percentOld!=percent){
        print("mse")
        sums <- rbind(sums, data.frame(percent = percentOld, sums = sumsVec))
        mse_mean <- rbind(mse_mean, data.frame(mseMean = mean((sumsVec-sum(bg_real$background))^2), 
                                               celltype=celltype, 
                                               percent = percentOld, 
                                               type=info[7], method=info[8]))
        data <- tmp
        percentOld = percent
        
        sumsVec <- c()
      } else {
        print("append")
        sumsVec <- c(sumsVec, sum(tmp$fraction))
        data <- rbind(data, tmp)
      }
      
    #}
  }
  mse_mean <- rbind(mse_mean, data.frame(mseMean = mean((sumsVec-sum(bg_real$background))^2), 
                                         celltype=celltype, 
                                         percent = percentOld, 
                                         type=info[7], method=info[8]))
  mse_mean <- rbind(mse_mean, data.frame(mseMean = mean(data$MSE), 
                                         celltype=celltype, 
                                         percent = percent, 
                                         type=info[7], method=info[8]))
  sums <- rbind(sums, data.frame(percent = percent, sums = sumsVec))
  return(list(background = bg_list, mse = mse_mean))
}
d = loadData(deconv_results_spikein)
mse_mean = d$mse
ggplot(mse_mean, aes(y=mseMean, x=percent))+geom_point()+geom_line()


g <- ggplot(rbind(d, data.frame(fraction = real$background, celltype=real$celltype, percent=0)))
scatter <- g+geom_point(aes(x=percent, y=fraction))
scatter
box <- g+geom_boxplot(aes(y=fraction, x=paste(percent)))
box
ggplot(d, aes(y=MSE, x=paste(percent)))+geom_boxplot()
ggplot(d, aes(y=meanMSE, x=percent))+geom_point()+geom_line()
