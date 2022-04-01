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
deconvolution_results <- list.files("/nfs/proj/omnideconv_benchmarking/pipelines/results_spikein/", pattern = "deconvolution_spikein.*Bx.xcell_.*tpm", full.names = TRUE)
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
  bg_dataframe <- NULL
  for(i in 1:length(deconv_results)){
    res <- deconv_results[i]
    print(res)
    result <- as.data.frame(readRDS(res))
    info <- unlist(strsplit(gsub(".rds", "", gsub("deconvolution_", "", basename(res))), "_"))
    celltype <- reEscapeCelltypes(info[4])
    fraction <- result[[celltype]]
    if(info[3]=="backgroundFor"){
      real <- data.frame(fraction = fraction, celltype=celltype,
                         type=info[5], method=info[6], spikein="background")
      #store as celltype_type_method
      bg_list[[paste(info[4], info[5], info[6], sep ="_")]] <- real
      bg_dataframe <- rbind(real, bg_dataframe)
    } else { #has 2 fields more
      percent <- as.integer(gsub("nrCells", "", info[5]))/1000
      tmp <- data.frame(fraction = fraction, celltype=celltype,
                        type=info[7], method=info[8], spikein = percent)
      #get by celltype_type_method
      tmp$background <- bg_list[[paste(info[4], info[7], info[8], sep ="_")]]$fraction
      data <- rbind(data, tmp)
    }
  }
  return(list(background_list = bg_list, background_df = bg_dataframe, predictions = data))#, mse = mse_mean))
}

d = loadData(deconvolution_results)

background_df = d$background_df
g <- ggplot(background_df, aes(fraction, celltype))+
  geom_boxplot()+
  facet_wrap("~method")+
  labs(x="background", y="cell type")
ggsave("spikeIn_background_noisePredictions_boxplot.jpeg", g)

all = d$predictions %>% 
  group_by(celltype, type, method, spikein) %>% 
  mutate(p_value=t.test(fraction, background, paired=FALSE, alternative = "greater")$p.value) %>%
  group_by(celltype, type, method, spikein, p_value) %>% 
  summarise_at(vars(fraction), 
               list(mean=mean, sd=sd, n=length))%>%
  mutate(ci=qt(0.975,df=n-1)*sd/sqrt(n)) %>% ungroup()

detection_limit = all %>% 
  group_by(method, celltype) %>% 
  #minimal fraction, at which the prediction is significantly different from the background. 
  filter(p_value < 0.05) %>% 
  summarise(minfrac = min(spikein))

g <- ggplot(all, aes(x=spikein, y=mean, color=method))+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+
  geom_point()+ facet_wrap("~celltype")+
  labs(x="fraction of spike-in cells", y="average estimate")
ggsave("spikeIn_avgEstimation.jpeg", g)

background_mean <- background_df %>% group_by(celltype, type, method, spikein) %>% 
  summarise_at(vars(fraction), list(mean_bg=mean))

g <- ggplot(all, aes(x=spikein, y=mean))+
  geom_ribbon(aes(ymin=mean-ci, ymax=mean+ci), alpha=.3) + 
  geom_point() + 
  facet_grid(method ~ celltype, scales = "free_y") + 
  geom_hline(mapping=aes(yintercept=0, color="zero"), linetype="dashed", show.legend=FALSE) + 
  geom_vline(mapping=aes(xintercept = 0, color="zero"), linetype="dashed", show.legend=FALSE) + 
  geom_hline(data = background_mean, mapping = aes(yintercept = mean_bg, color="false positive fraction")) +
  geom_vline(data = detection_limit, mapping = aes(xintercept = minfrac, color="detection limit")) +
  scale_color_manual(values=c("zero"="grey", "detection limit"="red", "false positive fraction"="blue"), 
                     guide=guide_legend("performance measure", override.aes = list(alpha=1)))+
  ylab("average estimate") + 
  xlab("fraction of spike-in cells")+ 
  theme(legend.position = "top")
ggsave("spikeIn_sensitivity_scatter.jpeg", g)

# 
# mse_mean = d$mse
# ggplot(mse_mean, aes(y=mseMean, x=percent))+geom_point()+geom_line()
# 
# 
# g <- ggplot(rbind(d, data.frame(fraction = real$background, celltype=real$celltype, percent=0)))
# scatter <- g+geom_point(aes(x=percent, y=fraction))
# scatter
# box <- g+geom_boxplot(aes(y=fraction, x=paste(percent)))
# box
# ggplot(d, aes(y=MSE, x=paste(percent)))+geom_boxplot()
# ggplot(d, aes(y=meanMSE, x=percent))+geom_point()+geom_line()
