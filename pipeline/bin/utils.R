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

# function to transform string to integer to use as seed
# from R package TeachingDemos by Greg Snow
char2seed <- function(x){

	tmp <- c(0:9,0:25,0:25)
	names(tmp) <- c(0:9,letters,LETTERS)

	x <- gsub("[^0-9a-zA-Z]","",as.character(x))

	xsplit <- tmp[ strsplit(x,'')[[1]] ]

	seed <- sum(rev( 7^(seq(along=xsplit)-1) ) * xsplit)
  seed <- as.integer( seed %% (2^31-1) )

	return(seed)
}

# Function to subset a dataset
subset_cells <- function(cell_matrix, annotations, batch_ids, num_cells, 
                         seed, coarse_annotations = NULL, fine_annotations = NULL){
  
  require(tidyverse)

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

# This function contains all the necessary steps to adapt omnideconv to the benchmarking workflow
signature_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, 
                                       method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path){
  
  # path to temporary directory that is inside results directory of one unique result 
  tmp_dir_path <- paste0(res_path, '/tmp/')

  # Signature building part
  # Dependent on which method we have 
  if (method == "cibersortx") {
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                          " 721a387e91c495174066462484674cb8") 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(tmp_dir_path,'_input')
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    cx_output <- paste0(tmp_dir_path,'_output')
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output)
    }

    signature <- omnideconv::build_model_cibersortx(
      sc_matrix, 
      annotations, 
      container = "docker", 
      verbose = TRUE, 
      input_dir = cx_input, 
      output_dir = cx_output
    )

  } else if (method == "dwls") {
    
    signature <- omnideconv::build_model_dwls(
      single_cell_object = sc_matrix,
      cell_type_annotations = annotations,
      dwls_method = "mast_optimized",
      path=res_path,
      diff_cutoff = 0.5,
      pval_cutoff = 0.05,
      verbose = TRUE,
      ncores = ncores
    )
    
  } else if (method == "scdc") {
    
    signature <- omnideconv::build_model_scdc(
      single_cell_object = sc_matrix,
      cell_type_annotations = annotations,
      batch_ids = sc_batch,
      markers = NULL,
      verbose = TRUE,
      ncores = ncores
    )$basis
    
  } else if(method == "scaden"){
    unlink(tmp_dir_path, recursive=TRUE)
    if(!dir.exists(paste0(tmp_dir_path))){
      dir.create(tmp_dir_path)
    }
    signature <- omnideconv::build_model_scaden(
      sc_matrix,
      annotations,
      bulk_matrix,
      temp_dir = tmp_dir_path,
      verbose = TRUE
    )

  }else if (method %in% c('autogenes', 'bayesprism', 'bisque', 'music')){
    signature <- NULL
      
  } else {
    message('Selected method is not supported in the benchmark. Please check again.')
    stop()
  }

  message('Finished signature building.')
  
  return(signature)
  
}

# This function contains all the necessary steps to adapt omnideconv to the benchmarking workflow
deconvolution_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, 
                                           signature, method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path,
                                           rmbatch_B_mode = FALSE, rmbatch_S_mode = FALSE){
  
  # path to temporary directory that is inside results directory of one unique result 
  tmp_dir_path <- paste0(res_path, '/tmp/')

  if(method=="autogenes"){

    deconvolution <- omnideconv::deconvolute_autogenes( 
      single_cell_object = sc_matrix,
      cell_type_annotations = annotations,
      bulk_gene_expression = bulk_matrix, 
      max_iter = 1000000,
      ngen = 5000,
      verbose = TRUE
    )$proportions
    colnames(deconvolution) <- reEscapeCelltypesAutogenes(colnames(deconvolution))

  } else if (method=="bayesprism"){
    
    # need to use different parameter in case of simulated bulk data
    if(grepl("simulation", bulk_name)){
      update_gibbs <- FALSE
      which_theta <- 'first'
    }else{
      update_gibbs <- TRUE
      which_theta <- 'final'
    }
    
    deconvolution <- omnideconv::deconvolute_bayesprism(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = annotations, 
      n_cores = ncores,
      species = args$species,
      update_gibbs = update_gibbs,
      which_theta = which_theta
    )$theta
    
  } else if (method=="bisque"){
    deconvolution <- t(omnideconv::deconvolute_bisque(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = annotations,
      batch_ids = sc_batch,
      verbose = TRUE,
    )$bulk.props)
    
  } else if (method=="cibersortx"){

    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                          " 721a387e91c495174066462484674cb8")  
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(tmp_dir_path,'_input')
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    cx_output <- paste0(tmp_dir_path,'_output')
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output)
    }

    # batch correction options for cibersortx
    #rmbatch_B_mode <- FALSE
    #rmbatch_S_mode <- TRUE
    #if(grepl("simulation" , bulk_name)){
    #  rmbatch_S_mode <- FALSE
    #}
    
    deconvolution <- omnideconv::deconvolute_cibersortx(
      single_cell_object = sc_matrix,
      bulk_gene_expression = bulk_matrix,
      cell_type_annotations = annotations, 
      signature = signature,
      container = 'docker',
      verbose = TRUE,
      input_dir = cx_input,
      output_dir = cx_output,
      rmbatch_B_mode = rmbatch_B_mode,
      rmbatch_S_mode = rmbatch_S_mode
    )
    unlink(cx_input, recursive=TRUE)
    unlink(cx_output, recursive=TRUE)
    
  } else if(method=="dwls") {
    
    deconvolution <- omnideconv::deconvolute_dwls(
        bulk_gene_expression = bulk_matrix, 
        signature = signature, 
        dwls_submethod = 'DampenedWLS',
        verbose = TRUE
      )

    #deconvolution <- NULL
    #tryCatch(    
    #  deconvolution <<- omnideconv::deconvolute_dwls(
    #    bulk_gene_expression = bulk_matrix, 
    #    signature = signature, 
    #    dwls_submethod = 'DampenedWLS',
    #    verbose = TRUE
    #  ), error = function(e){
    #    # DWLS can crash if there are not enough cells in the dataset or the cells are not differential enough between celltypes
    #    print(e[['message']])
    #    d <- matrix(
    #      rep(1, ncol(bulk_matrix) * ncol(signature)),
    #      nrow = ncol(bulk_matrix),
    #      ncol = ncol(signature)
    #    )
    #    deconvolution <<- data.frame(d)
    #    colnames(deconvolution) <<- as.character(colnames(signature))
    #    rownames(deconvolution) <<- as.character(colnames(bulk_matrix))
    #  }
    #)
    print(deconvolution)

  } else if (method == "music"){
    deconvolution <- omnideconv::deconvolute_music(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = annotations,
      batch_ids = sc_batch,
      verbose = TRUE,
    )$Est.prop.weighted

  } else if (method == "scdc") {    
    deconvolution <- omnideconv::deconvolute_scdc(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = annotations,
      batch_ids = sc_batch,
      verbose = TRUE,
    )$Est.prop.weighted

    if ("prop.est.mvw" %in% names(deconvolution)) {
      deconvolution <- deconvolution$prop.est.mvw
    } else if ("w_table" %in% names(deconvolution)) {
      deconvolution <- SCDC::wt_prop(deconvolution$w_table, deconvolution$prop.only)
    } else {
      message(
        "There seems to be an error, as the result of deconvolute_scdc did not ",
        "contain prop.est.mvw or w_table"
      )
    }
                
  } else if (method == "scaden") {
    if(!dir.exists(paste0(tmp_dir_path))){
      dir.create(tmp_dir_path)
    }
    deconvolution <- omnideconv::deconvolute_scaden(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      temp_dir = tmp_dir_path
    )
    unlink(tmp_dir_path, recursive=TRUE)

  } else {
    message('Selected method is not supported in the benchmark. Please check again.')
    stop()
  }
  
  return(deconvolution)
}

# This function and the excel file will be needed for the reannotation of the ground truth fractions
reannotate_facs <- function(facs.table, annotation, new.annotation.level){
  
  facs.table <- as.data.frame(facs.table)
  annotation <- annotation[which(annotation$Fine %in% colnames(facs.table)), ]
  cell_types <- unique(annotation[[new.annotation.level]])
  for(c in cell_types){
    # These are the fine cell types
    cur.cell.types <- annotation[which(annotation[[new.annotation.level]] == c), 1]
    
    if(length(cur.cell.types) > 1){
      facs.table[c] <- rowSums(facs.table[, c(cur.cell.types)])
      facs.table[, c(cur.cell.types)] <- NULL
    }
  }
  facs.table
}

reannotate_cell_types <- function(celltypes_fine, annotation, new.annotation.level){
  
  annotation <- annotation[which(annotation$Fine %in% celltypes_fine), ]
  cell_types <- unique(annotation[[new.annotation.level]])

  cell_types
}