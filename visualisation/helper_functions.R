get_performance_metrics <- function(output_dir, ct_values, query_sc, method_parameter_df){
  res <- lapply(ct_values, function(i) {
    l <- list.files(output_dir, pattern=paste0('*ct',i,'_'))
    
    ct_res <- lapply(l, function(f){
      if(grepl('vanderbilt_lung', f)){
        # need to replace the '_' in bulk ds name with '-' so that indexing is still correct
        f_new <- sub("^(([^_]*_){3}[^_]*)_([^_]*)", "\\1-\\3", f)
        xl <- strsplit(f_new, '_')[[1]]
      }else{
        xl <- strsplit(f, '_')[[1]]
      }
      method <- xl[1]
      sc_ds <- xl[2]
      sc_norm <- xl[3]
      ct <- as.character(gsub('ct','',xl[6]))
      
      if(method %in% method_parameter_df$method & sc_ds == query_sc & ct %in% ct_values){
        query_sc_norm <- method_parameter_df[which(method_parameter_df$method == method),'sc_norm']
        
        if(sc_norm == query_sc_norm){
          results_file <- paste0(output_dir,f,'/results_metric.rds')
          if(file.exists(results_file)){
            x <- readRDS(results_file)
            df <- x[['rmse_cell_type']]
            df <- cbind(df, x$cor_cell_type[,c(2,3)])
            df$method <- method
            df$sc_ds <- sc_ds
            df$sc_norm <- sc_norm
            df$bulk_ds <- xl[4]
            df$bulk_norm <- xl[5]
            df$ct <- ct
            df$replicate <- as.character(gsub('rep','',xl[7]))
            df$method_norm_combi <- paste0(method, sc_norm, df$bulk_norm)
            
            df
          }
        }
      }
    })
    
    df <- rbindlist(ct_res)
  })
  return(rbindlist(res))
}

get_performance_metrics_per_sample <- function(output_dir, ct_values, query_sc, method_parameter_df){
  res <- lapply(ct_values, function(i) {
    l <- list.files(output_dir, pattern=paste0('*ct',i,'_'))
    
    ct_res <- lapply(l, function(f){
      if(grepl('vanderbilt_lung', f)){
        # need to replace the '_' in bulk ds name with '-' so that indexing is still correct
        f_new <- sub("^(([^_]*_){3}[^_]*)_([^_]*)", "\\1-\\3", f)
        xl <- strsplit(f_new, '_')[[1]]
      }else{
        xl <- strsplit(f, '_')[[1]]
      }
      method <- xl[1]
      sc_ds <- xl[2]
      sc_norm <- xl[3]
      ct <- as.character(gsub('ct','',xl[6]))
      
      if(method %in% method_parameter_df$method & sc_ds == query_sc & ct %in% ct_values){
        query_sc_norm <- method_parameter_df[which(method_parameter_df$method == method),'sc_norm']
        
        if(sc_norm == query_sc_norm){
          results_file <- paste0(output_dir,f,'/results_metric.rds')
          if(file.exists(results_file)){
            x <- readRDS(results_file)
            df <- x[['rmse_sample']]
            df <- cbind(df, x$cor_sample[,c(2,3)])
            df$method <- method
            df$sc_ds <- sc_ds
            df$sc_norm <- sc_norm
            df$bulk_ds <- xl[4]
            df$bulk_norm <- xl[5]
            df$ct <- ct
            df$replicate <- as.character(gsub('rep','',xl[7]))
            df$method_norm_combi <- paste0(method, sc_norm, df$bulk_norm)
            
            df
          }
        }
      }
    })
    
    df <- rbindlist(ct_res)
  })
  return(rbindlist(res))
}

get_all_performance_metrics <- function(output_dir, ct_values, query_sc, method_parameter_df){
  res <- lapply(ct_values, function(i) {
    l <- list.files(output_dir, pattern=paste0('*ct',i,'_'))
    
    ct_res <- lapply(l, function(f){
      if(grepl('vanderbilt_lung', f)){
        # need to replace the '_' in bulk ds name with '-' so that indexing is still correct
        f_new <- sub("^(([^_]*_){3}[^_]*)_([^_]*)", "\\1-\\3", f)
        xl <- strsplit(f_new, '_')[[1]]
      }else{
        xl <- strsplit(f, '_')[[1]]
      }
      
      method <- xl[1]
      sc_ds <- xl[2]
      sc_norm <- xl[3]
      ct <- as.character(gsub('ct','',xl[6]))
      
      if(method %in% method_parameter_df$method & sc_ds == query_sc & ct %in% ct_values){
        query_sc_norm_vect <- as.vector(unlist(method_parameter_df[which(method_parameter_df$method == method),'sc_norm']))
        
        cur.dataframe <- data.frame()
        
        for(query_sc_norm in query_sc_norm_vect){
          if(sc_norm == query_sc_norm){
            results_file <- paste0(output_dir,f,'/results_metric.rds')
            if(file.exists(results_file)){
              x <- readRDS(results_file)
              df <- x[['rmse_cell_type']]
              df <- cbind(df, x$cor_cell_type[,c(2,3)])
              df$method <- method
              df$sc_ds <- sc_ds
              df$sc_norm <- sc_norm
              df$bulk_ds <- xl[4]
              df$bulk_norm <- xl[5]
              df$ct <- ct
              df$replicate <- as.character(gsub('rep','',xl[7]))
              df$method_norm_combi <- paste0(method, sc_norm, df$bulk_norm)
              
              cur.dataframe<-rbind(cur.dataframe, df)
            }
          }
        }
        
        cur.dataframe
      }
    })
    
    df <- rbindlist(ct_res)
  })
  return(rbindlist(res))
}

get_runtimes <- function(output_dir, sc_dir, preprocess_dir, ct_values, query_sc, methods){
  res <- lapply(ct_values, function(i) {
    l <- list.files(output_dir, pattern=paste0('*ct',i,'_'))
    
    ct_res <- lapply(l, function(f){
      xl <- strsplit(f, '_')[[1]]
      method <- xl[1]
      sc_ds <- xl[2]
      ct <- as.character(gsub('ct','',xl[6]))
      rep <- as.character(gsub('rep','',xl[7]))
      if(i == 0){
        sc_data_dir <- paste0(sc_dir, query_sc)
      }else{
        sc_data_dir <- paste0(preprocess_dir, query_sc,'_counts_perc',ct,'_rep',rep)
      }
      
      if(method %in% methods & sc_ds == query_sc & ct %in% ct_values){
        results_file <- paste0(output_dir,f,'/results_metric.rds')
        cell_anno <- readRDS(paste0(sc_data_dir,'/celltype_annotations.rds'))
        n_cells <- length(cell_anno)
        if(file.exists(results_file)){
          x <- readRDS(results_file)
          df <- x[['runtimes']]
          df$n_cells <- n_cells 
          df
        }
      }
    })
    
    df <- rbindlist(ct_res)
  })
  df <- rbindlist(res)
  df$ct_fraction <- as.character(df$ct_fraction)
  df
}

get_fractions <- function(output_dir, ct_values, query_sc, method_parameter_df){
  res <- lapply(ct_values, function(i) {
    l <- list.files(output_dir, pattern=paste0('*ct',i,'_'))
    
    ct_res <- lapply(l, function(f){
      if(grepl('vanderbilt_lung', f)){
        # need to replace the '_' in bulk ds name with '-' so that indexing is still correct
        f_new <- sub("^(([^_]*_){3}[^_]*)_([^_]*)", "\\1-\\3", f)
        xl <- strsplit(f_new, '_')[[1]]
      }else{
        xl <- strsplit(f, '_')[[1]]
      }
      method <- xl[1]
      sc_ds <- xl[2]
      sc_norm <- xl[3]
      ct <- as.character(gsub('ct','',xl[6]))
      if (method %in% method_parameter_df$method & sc_ds == query_sc & ct %in% ct_values){
        query_sc_norm <- method_parameter_df[which(method_parameter_df$method == method),'sc_norm']
        
        if(sc_norm == query_sc_norm){
          results_file <- paste0(output_dir,f,'/results_metric.rds')
          
          if(file.exists(results_file)){
            res <- readRDS(results_file)
            df_deconv <- melt(as.data.table(res$deconv.results, keep.rownames = T), id.vars = 1)
            colnames(df_deconv) <- c('sample','celltype','fraction')
            df_true <- melt(as.data.table(res[[5]], keep.rownames = T), id.vars = 1)
            colnames(df_true) <- c('celltype','sample','fraction')
            df <- merge(df_deconv, df_true, by=c('celltype','sample'), suffixes = c('.estimate','.true'))
            df$method <- method
            df$sc_ds <- sc_ds
            df$sc_norm <- sc_norm
            df$bulk_ds <- xl[4]
            df$bulk_norm <- xl[5]
            df$ct <- ct
            df$replicate <- as.character(gsub('rep','',xl[7]))
            df$method_norm_combi <- paste0(method, sc_norm, df$bulk_norm)
            
            df
          } 
        }
      }
    })
    
    df <- rbindlist(ct_res)
  })
  return(rbindlist(res))
}

get_all_fractions <- function(output_dir, ct_values, query_sc, method_parameter_df){
  res <- lapply(ct_values, function(i) {
    l <- list.files(output_dir, pattern=paste0('*ct',i,'_'))
    
    ct_res <- lapply(l, function(f){
      if(grepl('vanderbilt_lung', f)){
        # need to replace the '_' in bulk ds name with '-' so that indexing is still correct
        f_new <- sub("^(([^_]*_){3}[^_]*)_([^_]*)", "\\1-\\3", f)
        xl <- strsplit(f_new, '_')[[1]]
      }else{
        xl <- strsplit(f, '_')[[1]]
      }
      method <- xl[1]
      sc_ds <- xl[2]
      sc_norm <- xl[3]
      ct <- as.character(gsub('ct','',xl[6]))
      if (method %in% method_parameter_df$method & sc_ds == query_sc & ct %in% ct_values){
        query_sc_norm_vect <- as.vector(unlist(method_parameter_df[which(method_parameter_df$method == method),'sc_norm']))
        
        cur.dataframe <- data.frame()
        
        for(query_sc_norm in query_sc_norm_vect){
          
          
          if(sc_norm == query_sc_norm){
            results_file <- paste0(output_dir,f,'/results_metric.rds')
            
            if(file.exists(results_file)){
              res <- readRDS(results_file)
              df_deconv <- melt(as.data.table(res$deconv.results, keep.rownames = T), id.vars = 1)
              colnames(df_deconv) <- c('sample','celltype','fraction')
              df_true <- melt(as.data.table(res$facs_ground_truth, keep.rownames = T), id.vars = 1)
              colnames(df_true) <- c('celltype','sample','fraction')
              df <- merge(df_deconv, df_true, by=c('celltype','sample'), suffixes = c('.estimate','.true'))
              df$method <- method
              df$sc_ds <- sc_ds
              df$sc_norm <- sc_norm
              df$bulk_ds <- xl[4]
              df$bulk_norm <- xl[5]
              df$ct <- ct
              df$replicate <- as.character(gsub('rep','',xl[7]))
              df$method_norm_combi <- paste0(method, sc_norm, df$bulk_norm)
              
              cur.dataframe<-rbind(cur.dataframe, df)
            } 
          }
        }
        cur.dataframe
      }
    })
    
    df <- rbindlist(ct_res)
  })
  return(rbindlist(res))
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      iqr = quantile(x[[col]], probs=c(.25, .75), na.rm=TRUE))
  }
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- plyr::rename(data_sum, c("median" = varname))
  return(data_sum)
}

delta_correlation <- function(df){
  res <- unlist(lapply(1:nrow(df), function(i){
    ct_i <- which(levels(df$ct[i]) == df$ct[i])
    if(ct_i == 1){
      delta_correlation_i <- 0
    }else{
      ct_prev <- levels(df$ct)[1]
      # take 'current' correlation and subtract correlation at subset size with 1st position in levels of ct
      subset_df <- df %>% subset(cell_type == df$cell_type[i] &
                                   method == df$method[i] &
                                   bulk_ds == df$bulk_ds[i] &
                                   ct == ct_prev)
      if(nrow(subset_df) > 0){
        delta_correlation_i <- df$cor[i] - subset_df%>% pull(cor)
      }else{
        delta_correlation_i <- NA
      }
      
    } 
    return(delta_correlation_i)
  }))
}

delta_rmse <- function(df){
  res <- unlist(lapply(1:nrow(df), function(i){
    ct_i <- which(levels(df$ct[i]) == df$ct[i])
    if(ct_i == 1){
      delta_RMSE_i <- 0
    }else{
      ct_prev <- levels(df$ct)[1]
      # take 'current' RMSE and subtract RMSE at subset size with 1st position in levels of ct
      subset_df <- df %>% subset(cell_type == df$cell_type[i] &
                                   method == df$method[i] &
                                   bulk_ds == df$bulk_ds[i] &
                                   ct == ct_prev)
      if(nrow(subset_df) > 0){
        delta_RMSE_i <- df$RMSE[i] - subset_df%>% pull(RMSE)
      }else{
        delta_RMSE_i <- NA
      }
    }
    return(delta_RMSE_i)
  }))
}

plot_runtime <- function(runtimes_df, bulk_ds_name, log=FALSE){
  combined_runtimes <- runtimes_df[,sum(elapsed),by=c('method','replicate','ct_fraction','sc_ds','cs_norm','bulk_ds','bulk_norm','n_cells')]
  colnames(combined_runtimes)[which(colnames(combined_runtimes) == 'V1')] <- 'elapsed'
  combined_runtimes$process <- 'COMBINED'
  runtime_df <- rbind(runtimes_df, combined_runtimes, fill=TRUE)
  
  if(log){
    runtime_df$elapsed <- log(runtime_df$elapsed)
    runtime_df[which(elapsed == -Inf),'elapsed'] <- 0
    runtime_df$n_cells <- log(as.numeric(runtime_df$n_cells))
  }
  
  df <- data_summary(subset(runtime_df, bulk_ds == bulk_ds_name), 'elapsed', c('n_cells','process','method'))
  df$process <- factor(df$process, levels = c('SIGNATURE','DECONVOLUTION','COMBINED'))
  #df$n_cells <- factor(df$n_cells, levels = c(n_cells))
  
  ggplot(df, aes(x=n_cells, y=elapsed, group=method, color=method))+
    geom_point(alpha=.7)+
    geom_line(alpha=.7)+
    facet_wrap(~process)+
    geom_errorbar(aes(ymin=elapsed-sd, ymax=elapsed+sd), width=.2)+
    theme_bw()+
    ylab('log(elapsed time (s))')+
    xlab('log(number of single cells in reference)')+
    scale_color_brewer(palette = 'Set2')
}

plot_performance_replicates <- function(performance_df, bulk_ds_name, score, method_parameter_df){
  
  df <- data_summary(subset(performance_df, bulk_ds == bulk_ds_name & method_norm_combi %in% method_parameter_df$method_norm_combi), score, c('ct','cell_type','method'))
  colnames(df)[which(colnames(df) == score)] <- 'score'
  df$ct <- factor(df$ct, levels = ct_values)
  
  ggplot()+
    #geom_rect(data=df, aes(ymin=-Inf,ymax=0,xmin=-Inf,xmax=Inf), alpha=.01, fill='#CCD6DB')+
    #geom_rect(data=df, aes(ymin=0,ymax=0.5,xmin=-Inf,xmax=Inf), alpha=.01, fill='#e5eaed')+
    geom_point(data=df, aes(x=ct, y=score, group=method, color=method), alpha=.8)+
    geom_line(data=df, aes(x=ct, y=score, group=method, color=method), alpha=.8)+
    geom_errorbar(data=df, aes(x=ct, y=score, group=method, color=method, ymin=score-sd, ymax=score+sd), width=.2, alpha=.8)+
    facet_wrap(~cell_type, scales='free_x')+
    theme_bw()+
    ylab(paste0(score))+
    xlab('subset of cell type')+
    ggtitle(paste0(score, ' for different cell types, deconvoluted on ', bulk_ds_name, ' dataset; \nsignature created with increasing subset sizes of hao dataset'))+
    scale_color_brewer(palette = 'Set2')
}

plot_performance_replicates_scatter <- function(performance_df, bulk_ds_name, method_parameter_df){
  
  df_cor <- data_summary(subset(performance_df, bulk_ds == bulk_ds_name & method_norm_combi %in% method_parameter_df$method_norm_combi), 'cor', c('ct','cell_type','method'))
  colnames(df_cor)[which(colnames(df_cor) == 'sd')] <- 'sd.cor'
  df_rmse <- data_summary(subset(performance_df, bulk_ds == bulk_ds_name & method_norm_combi %in% method_parameter_df$method_norm_combi), 'RMSE', c('ct','cell_type','method'))
  colnames(df_rmse)[which(colnames(df_rmse) == 'sd')] <- 'sd.rmse'
  
  df <- cbind(df_cor, df_rmse[,c('RMSE','sd.rmse')])
  df$ct <- factor(df$ct, levels = ct_values)
  
  ggplot(df, aes(x=cor, y=RMSE, color=cell_type, group=cell_type))+
    geom_point()+
    #geom_line()+
    geom_errorbar(aes(ymin=RMSE-sd.rmse, ymax=RMSE+sd.rmse), width=.001)+
    geom_errorbar(aes(xmin=cor-sd.cor, xmax=cor+sd.cor), width=.001)+
    facet_grid(ct~method)+
    theme_bw()+
    scale_color_brewer(palette = 'Set2')
}

plot_performance_heatmap <- function(performance_df, bulk_ds_name, score, method_parameter_df){
  df <- data_summary(subset(performance_df, bulk_ds == bulk_ds_name & method_norm_combi %in% method_parameter_df$method_norm_combi), score, c('ct','cell_type','method'))
  colnames(df)[which(colnames(df) == score)] <- 'score'
  df$ct <- factor(df$ct, levels = ct_values)
  
  ggplot(df, aes(x=method, y=ct, fill=score))+
    geom_tile()+
    facet_wrap(~cell_type)+
    theme_bw()+scale_fill_gradient(low = 'white', high='#2d87bb')+
    geom_text(aes(label=round(score, digits=2), color=ifelse(score > 0.6, 'white','black')), size=2.5)+
    ylab('Number of single cells in reference')+
    scale_colour_manual(values=c("white"="white", "black"="black"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

plot_fractions <- function(fractions_df, bulk_name, sample_name=0){
  if(sample_name == 0){
    df <- data_summary(subset(fractions_df, bulk_ds == bulk_ds_name), 'value', c('ct','variable','method'))
  }else{
    df <- data_summary(subset(fractions_df, bulk_ds == bulk_ds_name & rn == sample_name), 'value', c('ct','variable','method'))
  }
  df$ct <- factor(df$ct, levels = ct_values)
  
  ggplot(df, aes(x=ct, y=value, group=method, color=method))+
    geom_point()+
    geom_line()+
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2)+
    facet_wrap(~variable, scales='free_x')+
    theme_bw()+
    ylab('estimated cell-type fraction')+
    xlab('subset of cell type')+
    ggtitle(paste0('Estimated cell-type fractions by methods over different subset sizes \n[for ',ifelse(sample_name==0,'all samples',paste0('single sample: ',sample_name)),']'))+
    scale_color_brewer(palette = 'Set2')
}

plot_scatter_replicates <- function(fractions_df, bulk_ds_name, method_parameter_df, m){
  df <- subset(fractions_df, method == m & bulk_ds == bulk_ds_name & method_norm_combi %in% method_parameter_df$method_norm_combi)
  df_summary <- data_summary(df, 'fraction.estimate', c('ct','celltype','sample','fraction.true'))
  df_summary$ct <- factor(df_summary$ct, levels = ct_values)
  
  ggplot(df_summary, aes(x=fraction.true, y=fraction.estimate))+
    geom_point(aes(color=celltype))+
    facet_wrap(~ct)+
    geom_errorbar(aes(group=sample, color=celltype, ymin=(fraction.estimate)-sd, ymax=(fraction.estimate)+sd), width=0, alpha=.8)+
    geom_abline(alpha=.75)+
    theme_bw()+
    ylab('Mean estimated cell-type fraction')+
    xlab('true cell-type fraction')+
    ggtitle(paste0('Comparison of true and estimated cell-type fractions with ', m,'\n for different subset sizes. Bulk dataset: ',bulk_ds_name))+
    scale_color_brewer(palette = 'Set2')+
    stat_cor(size=2.7)
}

plot_scatter <- function(fractions_df, bulk_ds_name, bulk_norm, m=NULL){
  if(is.null(m)){
    df <- subset(fractions_df, bulk_ds == bulk_ds_name & bulk_norm == bulk_norm)
  }else{
    df <- subset(fractions_df, bulk_ds == bulk_ds_name & bulk_norm == bulk_norm & method %in% m)
  }
  
  df$fraction.estimate[which(is.na(df$fraction.estimate))] <- 0
  df$fraction.true[which(is.na(df$fraction.true))] <- 0
  
  ggplot(df)+
    geom_point(aes(x=fraction.true, y=fraction.estimate, color=method))+
    facet_wrap(~celltype, scales='free')+
    geom_abline()+
    scale_color_brewer(palette = 'Set2')+
    theme_bw()
}

plot_distance_from_identity <- function(fractions_df, bulk_name, m){
  df <- subset(fractions_df, method == m & bulk_ds == bulk_name)
  df$distance_from_identity <- df$fraction.true - df$fraction.estimate
  df_summary <- data_summary(df, 'distance_from_identity', c('ct','celltype','sample'))
  df_summary$ct <- factor(df_summary$ct, levels = ct_values)
  
  
}
