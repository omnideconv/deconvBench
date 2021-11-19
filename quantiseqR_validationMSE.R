### compare results of quantiseq and quantiseqr ###
load("/nfs/data/omnideconv_benchmarking/quantiseqComparison/results_finotello.RData") #RData object in which all finotello results are stored
load("/nfs/data/omnideconv_benchmarking/quantiseqComparison/results_hoek.RData") #RData object in which all hoek results are stored
load("/nfs/data/omnideconv_benchmarking/data/PBMC/finotello/finotello_pbmc_facs.RData")#finotello facs
load("/nfs/data/omnideconv_benchmarking/data/PBMC/hoek/hoek_pbmc_facs.RData")#hoek facs

remapCelltypesHoek <- function(celltype){
  if(grepl(pattern = "NK", x = celltype)){
    return("NK.cells")
  } else if(grepl(pattern = "Bcells", x = celltype)){
    return("B.cells")
  } else if(grepl(pattern = "Tcell", x = celltype)){
    return("T.cells")
  } else if(grepl(pattern = "mono", x = celltype)){
    return("Monocytes")
  } else if(grepl(pattern = "DC", x = celltype)){
    return("Dendritic.cells")
  } else {
    return(celltype)
  }
}
### melt dfs and merge them into one ###
### returns df with rows of MSE, MSE per celltype, correlation and correlation for each celltype and the according value pairs###
calcCorStats <- function(arrays, tumor, scale, method, dataset){
  name <- paste(substr(arrays, 1, 1), substr(tumor, 1, 1), substr(scale, 1, 1), "_", method, "_", dataset, sep="")
  quantiseq <- get(paste("quantiseq", name, sep="_"))
  quantiseqr <- get(paste("quantiseqr", name, sep="_"))
  facs <- get(paste(dataset, "pbmc_facs", sep="_")) %>% rownames_to_column("celltype") %>% tidyr::pivot_longer(cols = !celltype, names_to = "Sample", values_to = "true_fraction")
  if(dataset=="hoek"){
    facs$celltype <- sapply(facs$celltype, function(x) remapCelltypesHoek(x))
    quantiseq$T.cells <- quantiseq$T.cells.CD4 + quantiseq$T.cells.CD8 + quantiseq$Tregs
    quantiseqr$T.cells <- quantiseqr$T.cells.CD4 + quantiseqr$T.cells.CD8 + quantiseqr$Tregs
  }
  tmp <- tidyr::pivot_longer(quantiseq, cols= !Sample, names_to = "celltype", values_to = "quantiseq_frac") %>% 
    full_join(pivot_longer(quantiseqr, cols= !Sample, names_to = "celltype", values_to = "quantiseqr_frac"), keys = c("Sample", "celltype")) %>% 
    left_join(facs, keys = c("Sample", "celltype"))%>% 
    group_by(celltype)
  returnDF <- data.frame(celltype=sort(unique(tmp$celltype)),
                     MSE = mean((tmp$quantiseq_frac-tmp$quantiseqr_frac)^2), 
                     MSE_celltypes = (tmp %>% group_by(celltype) %>% summarise(MSE_cell = mean((quantiseq_frac - quantiseqr_frac)^2)))$MSE_cell, 
                     correlation = cor(tmp$quantiseq_frac, tmp$quantiseqr_frac), 
                     correlation_celltypes = (tmp %>% group_by(celltype) %>% summarise(correlation_cell = cor(quantiseq_frac, quantiseqr_frac)))$correlation_cell, 
                     correlation_quantiseq_gtruth = (tmp %>% group_by(celltype) %>% summarise(correlation_cell = cor(quantiseq_frac, true_fraction)))$correlation_cell, 
                     correlation_quantiseqr_gtruth = (tmp %>% group_by(celltype) %>% summarise(correlation_cell = cor(quantiseqr_frac, true_fraction)))$correlation_cell) %>%
    mutate(arrays=arrays, tumor=tumor, scale=scale, method=method, .before=celltype)
  return(returnDF)
}

### returns a table with MSE for all param combinations ###
makeMSEtable <- function(dataset, celltype_wise=TRUE){
  final <- crossing(arrays = c(TRUE, FALSE), tumor = c(TRUE, FALSE), scale = c(TRUE, FALSE), method = c("lsei", "hampel", "huber", "bisquare"))
  full <- NULL
  for(i in 1:(nrow(final))){
      full = rbind(full, calcCorStats(final[i,]$arrays, final[i,]$tumor, final[i,]$scale, final[i,]$method, dataset))
    }
  return(full) 
}


# retrieve dataframes and save them 
hoekMSEtable <- makeMSEtable("hoek")
write_tsv(hoekMSEtable, "/nfs/data/omnideconv_benchmarking/quantiseqComparison/hoekMSEtable.tsv")
finotelloMSEtable <- makeMSEtable("finotello")
write_tsv(finotelloMSEtable, "/nfs/data/omnideconv_benchmarking/quantiseqComparison/finotelloMSEtable.tsv")
