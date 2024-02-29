library(Seurat)
library(tidyverse)

escapeCelltypesAutogenes <- function(celltype){
  celltype <- gsub("\\+", "21b2c6e87f8711ec9bf265fb9bf6ab9c", celltype)
  celltype <- gsub("-", "21b2c7567f8711ec9bf265fb9bf6ab9a", celltype)
  celltype <- gsub("\\(", "21b2c7567f8711ec9bf265fb9bf6ab9f", celltype)
  celltype <- gsub(")", "21b2c7567f8711ec9bf265fb9bf6ab9g", celltype)
  return(gsub(" ", "xxxx", celltype))
}
reEscapeCelltypesAutogenes <- function(celltype){ 
  return(gsub("xxxx", " ", celltype))
}


# Function to subset a dataset
subset_cells <- function(cell_matrix, annotations, batch_ids, num_cells, seed, coarse_annotations = NULL, fine_annotations = NULL){
  
  seurat.obj <- Seurat::CreateSeuratObject(counts=cell_matrix,
                                           assay="RNA")
  seurat.obj <- Seurat::AddMetaData(seurat.obj, 
                            batch_ids, 
                            'batch_ids')
  seurat.obj <- Seurat::AddMetaData(seurat.obj,
                            annotations,
                            'cell_type')
  set.seed(seed)

  if(is.null(coarse_annotations)){
    sampled.metadata <- seurat.obj@meta.data %>%
      rownames_to_column(., 'barcode') %>%
      group_by(., cell_type) %>% 
      nest() %>%            
      mutate(n =  map_dbl(data, nrow)) %>%
      mutate(n = min(n, num_cells)) %>%
      ungroup() %>% 
      mutate(samp = map2(data, n, sample_n)) %>% 
      select(-data) %>%
      unnest(samp)

      sampled.data <- subset(seurat.obj, 
                             cells = sampled.metadata$barcode)
      ls <- list('data' = sampled.data[["RNA"]]$counts %>%
                  as.matrix(.),
                 'annotations' = sampled.data@meta.data$cell_type,
                 'batch_id' = sampled.data@meta.data$batch_ids)
  
  } else {
      seurat.obj <- Seurat::AddMetaData(seurat.obj,
                                coarse_annotations,
                                'cell_type_coarse')
      seurat.obj <- Seurat::AddMetaData(seurat.obj,
                                fine_annotations,
                                'cell_type_fine')   
      sampled.metadata <- seurat.obj@meta.data %>%
        rownames_to_column(., 'barcode') %>%
        group_by(., cell_type_fine) %>% 
        nest() %>%            
        mutate(n =  map_dbl(data, nrow)) %>%
        mutate(n = min(n, num_cells)) %>%
        ungroup() %>% 
        mutate(samp = map2(data, n, sample_n)) %>% 
        select(-data) %>%
        unnest(samp)     


        sampled.data <- subset(seurat.obj, 
                             cells = sampled.metadata$barcode)
        ls <- list('data' = sampled.data[["RNA"]]$counts %>%
                    as.matrix(.),
                   'annotations' = sampled.data@meta.data$cell_type,
                   'annotations_fine' = sampled.data@meta.data$cell_type_fine,
                   'annotations_coarse' = sampled.data@meta.data$cell_type_coarse,
                   'batch_id' = sampled.data@meta.data$batch_ids)                                            
  }

  ls
  
}
# This function contains all the necessary steps to adapt omnideconv to the benchmarking cluster
signature_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path, spillover = FALSE){
  
  # Signature building part
  # Dependent on which method we have 
  if (method %in% c('autogenes','bisque','music','scdc','bayesprism')){
    signature <- NULL
  } else if(method == 'scaden'){
    annotations <- escapeCelltypesAutogenes(annotations)
    
    signature_dir <- paste0(res_path, '/tmp/', method, '_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot/")
    
    unlink(signature_dir, recursive=TRUE)
    if(!dir.exists(paste0(signature_dir))){
      dir.create(signature_dir, recursive = TRUE)
    }
    
    scaden_tmp <- paste0(res_path, '/tmp/scaden_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    signature <- omnideconv::build_model_scaden(
        sc_matrix,
        annotations,
        bulk_matrix,
        temp_dir = scaden_tmp,
        verbose = TRUE
    )
  } else if(method == 'cibersortx'){
    
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                           " 721a387e91c495174066462484674cb8") 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(res_path, '/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input, recursive = TRUE)
    }
    
    cx_output <- paste0(res_path, '/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output, recursive = TRUE)
    }
    
    signature <- omnideconv::build_model_cibersortx(
      sc_matrix, 
      annotations, 
      container = "docker", 
      verbose = TRUE, 
      input_dir = cx_input, 
      output_dir = cx_output
    )
    
  } else if(method == 'dwls'){
    
    signature <- omnideconv::build_model_dwls(
      single_cell_object = sc_matrix,
      cell_type_annotations = annotations,
      dwls_method = "mast_optimized",
      path=res_path,
      diff_cutoff = 0.5, pval_cutoff = 0.05,
      verbose = TRUE,
      ncores = ncores
    )
    
  } else if (method == "rectangle") {
    AnnData <- reticulate::import("anndata")
    pd <- reticulate::import("pandas")
    # convert sc_matrix to adata object
    counts <- as.data.frame(t(sc_matrix))
    annotations_df <- pd$DataFrame(list(annotations))
    annotations_df <- t(annotations_df)
    rownames(annotations_df) <- rownames(counts)
    colnames(annotations_df) <- c("cell_type")
    annotations_df <- as.data.frame(annotations_df)
    
    adata <- AnnData$AnnData(counts, obs=annotations_df)
    row_indices <- rownames(bulk_matrix)

    # remove counts and annotations_df to free memory
    rm(counts, annotations_df)
    signature <- rp$pp$build_rectangle_signatures(adata,bulk_genes=row_indices)
    rectangle_tmp <- paste0('/vol/omnideconv_results/results_tmp/rectangle_',sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", annotation_category, "_annot")
    if(!dir.exists(rectangle_tmp)) {
      dir.create(rectangle_tmp)
    }
    # save signature to file as pkl
    reticulate::py_save_object(signature, file = paste0(rectangle_tmp, "/signature_result.pkl"))
  } else{
    signature <- omnideconv::build_model(
      single_cell_object = sc_matrix,
      cell_type_annotations = annotations,
      bulk_gene_expression = bulk_matrix,
      markers = NULL,
      batch_ids = sc_batch,
      verbose = TRUE,
      method = method
    )
  }
  
  #saveRDS(signature, file=paste0(res_path, "/signature.rds"))
  
  signature
  
}


deconvolution_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, signature, method, bulk_matrix, bulk_name, bulk_norm, rmbatch_B_mode = FALSE, rmbatch_S_mode = FALSE, ncores, res_path){
  
  if(method == 'autogenes'){
    deconvolution <- omnideconv::deconvolute_autogenes(
      single_cell_object = sc_matrix,
      cell_type_annotations = annotations, 
      bulk_gene_expression = bulk_matrix, 
      max_iter = 10000,
      ngen=5000,
      verbose=TRUE
    )$proportions
    colnames(deconvolution) <- reEscapeCelltypesAutogenes(colnames(deconvolution))
  } else if(method == 'scaden'){
    
    scaden_tmp <- paste0(res_path, '/tmp/scaden_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    deconvolution <- omnideconv::deconvolute_scaden(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      temp_dir = scaden_tmp
    )
    colnames(deconvolution) <- reEscapeCelltypesAutogenes(colnames(deconvolution))
    #unlink(scaden_tmp, recursive=TRUE)
    
  } else if(method == 'cibersortx'){
    
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                           " 721a387e91c495174066462484674cb8") 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(res_path, '/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input, recursive = TRUE)
    }
    
    cx_output <- paste0(res_path, '/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output, recursive = TRUE)
    }
    
    deconvolution <- omnideconv::deconvolute_cibersortx(
      bulk_gene_expression = bulk_matrix, 
      signature = signature, 
      container = 'docker',
      verbose = TRUE,
      input_dir = cx_input,
      output_dir = cx_output,
      rmbatch_B_mode = rmbatch_B_mode,
      rmbatch_S_mode = rmbatch_S_mode
    )
    #unlink(cx_input, recursive=TRUE)
    
  } else if(method == 'dwls'){
    
    #deconvolution <- NULL
    tryCatch(    
      deconvolution <<- omnideconv::deconvolute_dwls(
        bulk_gene_expression = bulk_matrix, 
        signature = signature, 
        dwls_submethod = 'DampenedWLS',
        verbose = TRUE
      ), error = function(e){
        print(e[['message']])
        d <- matrix(
          rep(1, ncol(bulk_matrix) * ncol(signature)),
          nrow = ncol(bulk_matrix),
          ncol = ncol(signature)
        )
        deconvolution <<- data.frame(d)
        #colnames(deconvolution) <<- as.character(colnames(signature))
        #rownames(deconvolution) <<- as.character(colnames(bulk_matrix))
      }
    )
    
  } else if(method == 'bayesprism'){
    
    deconvolution <- omnideconv::deconvolute_bayesprism(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = annotations, 
      n_cores = ncores,
      species = 'hs',
      update_gibbs = FALSE,
      which_theta = 'first'
    )$theta
    
  } else if (method == 'rectangle') {
    pd <- reticulate::import("pandas")
    AnnData <- reticulate::import("anndata")
    rectangle_tmp <- paste0('/vol/omnideconv_results/results_tmp/rectangle_',sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_ct", annotation_category, "_annot")
    if(!dir.exists(rectangle_tmp)) {
      dir.create(rectangle_tmp)
    }
    signature <- reticulate::py_load_object(paste0(rectangle_tmp, "/signature_result.pkl"))
    bulk_obs <- colnames(bulk_matrix)
    bulk_obs_df <- pd$DataFrame(list(bulk_obs))
    bulk_obs_df <- t(bulk_obs_df)
    rownames(bulk_obs_df) <- colnames(bulk_matrix)
    colnames(bulk_obs_df) <- c("bulk")
    bulk_obs_df <- as.data.frame(bulk_obs_df)
    bulk_matrix <- t(bulk_matrix)
    bulk_matrix <- as.data.frame(bulk_matrix)


    bulk_adata <- AnnData$AnnData(bulk_matrix, obs=bulk_obs_df)
    deconvolution <- rp$tl$deconvolution(signatures = signature, bulks = bulk_adata)

  } else {
    
    deconvolution <- omnideconv::deconvolute(
      bulk_gene_expression = bulk_matrix, 
      signature = signature, 
      single_cell_object = sc_matrix, 
      batch_ids = sc_batch, 
      cell_type_annotations = annotations, 
      method = method
    )
    
  }
  
  #colnames(deconvolution) <- gsub('\\.',' ',colnames(deconvolution))
  deconvolution
}

compute_metrics <- function(ref_facs, deconvolution, metric = c('cor', 'rmse'),
                            comparison = c('sample', 'cell_type')){
  
  compute_rmse <- function(x, y, zscored = FALSE, weighted = FALSE){
    if(zscored == TRUE){
      sd_x = sd(x)
      sd_y = sd(y)
      if(sd_x == 0){sd_x = 1}
      if(sd_y == 0){sd_y = 1}
      x <- (x - mean(x))/(sd_x)
      y <- (y - mean(y))/(sd_y)
    }
    
    rmse <- sqrt(mean((x - y)^2)) 
    if(weighted == TRUE){rmse <- rmse/max(y)}
    return(rmse)
  }
  
  
  
  deconvolution <- as.data.frame(deconvolution) %>%
    rownames_to_column(., 'sample') %>%
    gather(., key = 'cell_type', value = 'estimated_frac', -sample)
  
  ref_facs <- as.data.frame(ref_facs) %>%
    rownames_to_column(., 'cell_type') %>%
    gather(., key = 'sample', value = 'true_frac', -cell_type)
  
  results_df <- deconvolution %>%
    inner_join(., ref_facs) #, by = c('cell_type' ,'sample'))
  
  results_df$estimated_frac[which(is.na(results_df$estimated_frac))] <- 0
  results_df$true_frac[which(is.na(results_df$true_frac))] <- 0
  
  if(metric == 'cor'){
    res_metric <- results_df %>%
      group_by(!!! syms(comparison)) %>%
      dplyr::summarize(cor = cor.test(estimated_frac, true_frac)$estimate,
                       pval = cor.test(estimated_frac, true_frac)$p.value)
    
    cor_all_elements <- data.frame(comp = 'all' ,
                                   'cor' = cor.test(results_df$estimated_frac, results_df$true_frac)$estimate,
                                   'pval' = cor.test(results_df$estimated_frac, results_df$true_frac)$p.value,
                                   row.names = (nrow(res_metric) +1))
    colnames(cor_all_elements)[1] <- comparison
    res_metric <- rbind(res_metric, cor_all_elements)
    
  } else if(metric == 'rmse'){
    res_metric <- results_df %>%
      group_by(!!! syms(comparison)) %>%
      dplyr::summarize(RMSE = compute_rmse(estimated_frac, true_frac),
                       zRMSE = compute_rmse(estimated_frac, true_frac, zscored=TRUE), 
                       wRMSE = compute_rmse(estimated_frac, true_frac, weighted=TRUE)
                       )
    
    rmse_all_elements <- data.frame(comp = 'all' ,
                                    'RMSE' = compute_rmse(results_df$estimated_frac, results_df$true_frac),
                                    'zRMSE' = compute_rmse(results_df$estimated_frac, results_df$true_frac, zscored=TRUE),
                                    'wRMSE' = compute_rmse(results_df$estimated_frac, results_df$true_frac, weighted=TRUE),
                                    row.names = (nrow(res_metric) +1))
    
    colnames(rmse_all_elements)[1] <- comparison
    res_metric <- rbind(res_metric, rmse_all_elements)
  }
  
  return(res_metric)
  
}