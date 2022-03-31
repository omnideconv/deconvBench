#general plotting functions
library(forcats)
library(ggplot2)
library(dplyr)

#split up by one condition (eg sc set), shape by another (eg bulk), facet by method and celltype, color by celltype
separateByCondition <- function(filepath, data, subsetVectorFull, basename, shape, color="celltype", scales = "fixed", alpha = 0.2, addLM = FALSE){
  for(condition in unique(subsetVectorFull)){
    sub <- data[subsetVectorFull==condition,]
    corData <- drop_na(sub) %>% select(-sample) %>% group_by(facet, method) %>%  
      mutate(cor = round(cor(true_value, predicted_value), digits=2)) %>%
      mutate(rfill = ifelse(cor > r_threshold, higher_label, lower_label))
    corDataFill <- subset(corData, rfill==higher_label) %>% 
      distinct(facet, method, .keep_all = TRUE)
    corDataLabels <- corData %>% 
      distinct(facet, cor, method, .keep_all = TRUE) %>%
      mutate(x = 0.03, y = 0.4)
    p <- ggplot(corData, aes(x=true_value, y=predicted_value))+
      geom_rect(data = corDataFill, aes(fill = rfill),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,
                alpha = alpha, fill="yellow") +
      geom_point(aes(color=get(color), shape = get(shape)))+geom_abline()+
      facet_grid("method~facet", scales = scales)+
      scale_color_discrete(na.translate=F)+
      theme(axis.text.x = element_text(angle=90), legend.position = "bottom")+
      ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+
      labs(title = condition, shape = shape, color=color)
    if(addLM){
      p <- p + stat_smooth(method = "lm")
      p
      ggsave(filename = file.path(filepath, 
                                  paste(basename, "_scatter_lm_",  condition, ".jpeg", sep="")),
             p, width = 13, height=9)
    } else {
      p
      ggsave(filename = file.path(filepath, 
                                  paste(basename, "_scatter_",  condition, ".jpeg", sep="")),
             p, width = 13, height=9)
    }
    g <- ggplot(corDataLabels, aes(facet, fct_rev(method)))+
      geom_tile(aes(fill=cor))+
      geom_text(aes(label=round(cor, 2)))+
      scale_color_discrete(na.translate=F)+
      scale_fill_gradient(low = "white", high = "#1b98e0")+
      theme(axis.text.x = element_text(angle=90))+ 
      labs(title = condition, x="celltype", y="method")
    g
    ggsave(filename = file.path(filepath, 
                                paste(basename, "_correlation_", condition, ".jpeg", sep="")),
           g, width = 13, height=9)
  }
}


allTogether <- function(filepath, data, groupingVar1, groupingVar2, basename, shape, color="celltype", alpha = 0.2, addLM = FALSE, width = 13, height=9){
  corData <- drop_na(data) %>% select(-sample) %>% group_by(get(groupingVar1), get(groupingVar2)) %>%  
    mutate(cor = round(cor(true_value, predicted_value), digits=2)) %>%
    mutate(rfill = ifelse(cor > r_threshold, higher_label, lower_label))
  corDataFill <- subset(corData, rfill==higher_label) %>% 
    distinct(get(groupingVar1), get(groupingVar2), .keep_all = TRUE)
  corDataLabels <- corData %>% 
    distinct(get(groupingVar1), cor, get(groupingVar2), .keep_all = TRUE) %>%
    mutate(x = 0.03, y = 0.6)
  p <- ggplot(corData, aes(x=true_value, y=predicted_value))+
    geom_rect(data = corDataFill, aes(fill = rfill),
              xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,
              alpha = alpha, fill="yellow") +
    geom_point(aes(color=get(color), shape=get(shape)))+geom_abline()+
    ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+
    facet_grid(paste(groupingVar1, groupingVar2, sep="~"), labeller = label_wrap_gen(width=10))+
    scale_color_discrete(na.translate=F)+
    theme(axis.text.x = element_text(angle=90), legend.position = "bottom")+
    labs(shape=shape, color=color, X="true value", y="predicted value")
  if(addLM){
    p <- p + stat_smooth(method = "lm")
    p
    ggsave(filename = file.path(filepath, paste(basename, "_scatter_lm_group", groupingVar1, groupingVar2, ".jpeg", sep="")),p, width = width, height=height)
  } else {
    p
    ggsave(filename = file.path(filepath, paste(basename, "_scatter_group", groupingVar1, groupingVar2, ".jpeg", sep="")),p, width = width, height=height)
  }
  
  g <- ggplot(corDataLabels, aes(y=fct_rev(get(groupingVar1)), x=get(groupingVar2)))+
    geom_tile(aes(fill=cor))+
    geom_text(aes(label=round(cor, 2)))+
    scale_fill_gradient(low = "white", high = "#1b98e0")+
    scale_color_discrete(na.translate=F)+
    #scale_y_discrete(groupingVar1, trans="reverse")+
    theme(axis.text.x = element_text(angle=90)) + 
    labs(x=groupingVar2, y=groupingVar1)
  g
  ggsave(file.path(filepath, paste(basename, "_correlation_group", groupingVar1, groupingVar2, ".jpeg", sep="")), g, width = 13, height=9)
  
  subset <- data %>% group_by(combination, celltype) %>% 
    mutate(rmse = sqrt(mean((true_value - predicted_value)^2))) 
  #boxplot of RMSE with facets
  g <- ggplot(subset, aes(rmse, fct_rev(get(groupingVar1))))+
    geom_boxplot()+
    labs(x="RMSE", y = groupingVar1)+
    facet_wrap(~ get(groupingVar2))
  g
  ggsave(filename = file.path(filepath, paste(basename, "_boxplot_group", groupingVar1, groupingVar2, ".jpeg", sep="")),g, width = 10, height=7)
  
}

allTogetherOneFacetNoShape <- function(filepath, data, groupingVar1, basename, color="celltype", alpha = 0.2, addLM = FALSE){
  corData <- drop_na(data) %>% select(-sample) %>% group_by(get(groupingVar1)) %>%  
    mutate(cor = round(cor(true_value, predicted_value), digits=2)) %>%
    mutate(rfill = ifelse(cor > r_threshold, higher_label, lower_label))
  corDataFill <- subset(corData, rfill==higher_label) %>% 
    distinct(get(groupingVar1), .keep_all = TRUE)
  corDataLabels <- corData %>% 
    distinct(get(groupingVar1), cor, .keep_all = TRUE) %>%
    mutate(x = 0.03, y = 0.6)
  p <- ggplot(corData, aes(x=true_value, y=predicted_value))+
    geom_rect(data = corDataFill, aes(fill = rfill),
              xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,
              alpha = alpha, fill="yellow") +
    geom_point(aes(color=get(color)))+geom_abline()+
    ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=3, color="black", label.x.npc = 0.005)+
    facet_grid(~ get(groupingVar1))+
    scale_color_discrete(na.translate=F)+
    theme(axis.text.x = element_text(angle=90), legend.position = "bottom")+
    labs(color=color, X="true value", y="predicted value")
  if(addLM){
    p <- p + stat_smooth(method = "lm")
    p
    ggsave(filename = file.path(filepath, paste(basename, "_scatter_lm_group", groupingVar1, ".jpeg", sep="")),p, width = 16, height=7)
  } else {
    p
    ggsave(filename = file.path(filepath, paste(basename, "_scatter_group", groupingVar1, ".jpeg", sep="")),p, width = 16, height=7)
  }
  
  g <- ggplot(corDataLabels, aes(x=get(groupingVar1), y=y))+
    geom_tile(aes(fill=cor))+
    geom_text(aes(label=round(cor, 2)))+
    scale_fill_gradient(low = "white", high = "#1b98e0")+
    scale_color_discrete(na.translate=F)+
    #scale_y_discrete(groupingVar1, trans="reverse")+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    labs(y="", x=groupingVar1)
  g
  ggsave(file.path(filepath, paste(basename, "_correlation_group", groupingVar1, ".jpeg", sep="")), g, width = 16, height=7)
  
  subset <- data %>% group_by(combination, celltype) %>% 
    mutate(rmse = sqrt(mean((true_value - predicted_value)^2))) 
  #boxplot of RMSE with facets
  g <- ggplot(subset, aes(y=rmse, x=get(groupingVar1)))+
    geom_boxplot()+
    labs(y="RMSE", x = groupingVar1)
  g
  ggsave(filename = file.path(filepath, paste(basename, "_boxplot_group", groupingVar1, ".jpeg", sep="")),g, width = 10, height=7)
  
}

compareGroundTruthAllCombinations <- function(filepath, data, basename, scenario = "regular", shape = "sctype", color = "celltype", alpha = 0.2, addLM = FALSE){
  corData <- drop_na(data) %>% select(-sample) %>% group_by(combination, facet, celltype) %>%  
    mutate(cor = round(cor(true_value, predicted_value), digits=2)) %>%
    mutate(rfill = ifelse(cor > r_threshold, higher_label, lower_label))
  corDataLabels <- corData %>% 
    distinct(combination, facet, celltype, .keep_all = TRUE) %>%
    mutate(x = 0.03, y = 0.6)
  if(scenario=="regular"){
    corDataFill <- subset(corData, rfill==higher_label) %>% 
      distinct(facet, scset, bulk, .keep_all = TRUE)
    p <- ggplot(corData, aes(x=true_value, y=predicted_value))+
      geom_rect(data = corDataFill, aes(fill = rfill),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,
                alpha = alpha, fill="yellow") +
      geom_point(aes(color=get(color), shape=get(shape)))+geom_abline()+
      scale_color_discrete(na.translate=F)+
      ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+ 
      facet_wrap(facet ~ bulk + scset, ncol = length(unique(data$bulk))+length(unique(data$scset)))+
      theme(axis.text.x = element_text(angle=90), legend.position = "bottom")+
      labs(shape=shape, color=color)
  } else {
    corDataFill <- subset(corData, rfill==higher_label) %>% 
      distinct(facet, scset, .keep_all = TRUE)
    p <- ggplot(corData, aes(x=true_value, y=predicted_value))+
      geom_rect(data = corDataFill, aes(fill = rfill),
                xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,
                alpha = alpha, fill="yellow") +
      geom_point(aes(color=get(color), shape=get(shape)))+geom_abline()+
      scale_color_discrete(na.translate=F)+
      ggpubr::stat_cor(aes(label = ..r.label..), label.sep = "\n", size=2.5, color="black", label.x.npc = 0.005)+ 
      facet_wrap(facet ~ scset, ncol = length(unique(data$scset)))+
      theme(axis.text.x = element_text(angle=90), legend.position = "bottom")+
      labs(shape=shape, color=color)
  }
  if(addLM){
    p <- p + stat_smooth(method = "lm")
    p
    ggsave(filename = file.path(filepath, paste(basename, "_scatter_lm_allCombinations.jpeg", sep="")),
           p, width = 13, height=16)
  } else {
    p
    ggsave(filename = file.path(filepath, paste(basename, "_scatter_allCombinations.jpeg", sep="")),
           p, width = 13, height=16)
  }
  
  g <- ggplot(corDataLabels, aes(y=fct_rev(facet), x=combination))+
    geom_tile(aes(fill=cor))+
    geom_text(aes(label=round(cor, 2)))+
    scale_fill_gradient(low = "white", high = "#1b98e0")+
    scale_color_discrete(na.translate=F)+
    #scale_y_discrete(groupingVar1, trans="reverse")+
    theme(axis.text.x = element_text(angle=90)) + 
    labs(y="celltype", x="combination")
  g
  ggsave(file.path(filepath, paste(basename, "_correlation_allCombinations.jpeg", sep="")), 
         g, width = 13, height=9)
}
