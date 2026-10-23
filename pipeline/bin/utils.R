# normalizations that a method expects for the single-cell reference and the bulk data
method_normalization <- function(baseDir, method){
  normalizations <- read.table(paste0(baseDir, '/optimal_normalizations.csv'), sep = ',', header = TRUE)
  if(!method %in% normalizations$method){
    stop(paste0('No entry in optimal_normalizations.csv for method: ', method))
  }
  as.list(normalizations[normalizations$method == method, c('sc_norm', 'bulk_norm')])
}

# expression matrix of a single-cell dataset stored in <sc_path>/<sc_dataset>/
read_sc_matrix <- function(sc_path, sc_dataset, sc_norm){
  matrix_file <- if(sc_norm == 'counts') 'matrix_counts.rds' else 'matrix_norm_counts.rds'
  readRDS(file.path(sc_path, sc_dataset, matrix_file))
}

# expression matrix <bulk_dir>/<bulk_name>_<bulk_norm>.rds
read_bulk_matrix <- function(bulk_dir, bulk_name, bulk_norm){
  as.matrix(readRDS(file.path(bulk_dir, paste0(bulk_name, '_', bulk_norm, '.rds'))))
}

# elements of a list parameter that Nextflow passes as text, e.g. '[B cells, T cells CD8]'
parse_list_argument <- function(x, split = ',', trim = FALSE){
  elements <- strsplit(gsub('\\[|]', '', x), split)[[1]]
  if(trim){
    elements <- trimws(elements)
  }
  elements
}

# rectangle takes the single-cell reference as AnnData object
rectangle_reference <- function(sc_matrix, annotations){
  anndata::AnnData(X = t(sc_matrix),
                   obs = data.frame('cell_type' = annotations, row.names = colnames(sc_matrix)),
                   var = data.frame('genes' = rownames(sc_matrix), row.names = rownames(sc_matrix)))
}

# Settings of a method that the workflow passes as text, e.g. 'marker_fraction=0.1;other=x',
# as a named list. An absent argument gives an empty list.
parse_method_settings <- function(x){
  if(is.null(x) || is.na(x) || x == ''){
    return(list())
  }
  pairs <- strsplit(strsplit(x, ';')[[1]], '=')
  setNames(lapply(pairs, function(pair) pair[2]), sapply(pairs, function(pair) pair[1]))
}

# Groups of cell types for the hierarchical MuSiC (method 'music-hc'): a named vector
# cell type -> group, derived from two annotations of the same cells. Returns NULL if the coarser
# annotation does not group the cell types (a cell type in two groups, or no group with more than
# one cell type).
celltype_group_map <- function(annotations, coarser_annotations){
  if(length(annotations) != length(coarser_annotations)){
    return(NULL)
  }
  pairs <- unique(data.frame(cell_type = as.character(annotations), group = as.character(coarser_annotations), stringsAsFactors = FALSE))
  if(any(duplicated(pairs$cell_type)) || length(unique(pairs$group)) == nrow(pairs)){
    return(NULL)
  }
  setNames(pairs$group, pairs$cell_type)
}

# Group map for the annotation of a single-cell dataset, taken from the coarsest annotation level
# stored with it (celltype_annotations_coarse.rds, then celltype_annotations_normal.rds).
# Only read for 'music-hc'; NULL for every other method.
read_celltype_groups <- function(method, sc_dir, annotations){
  if(method != 'music-hc'){
    return(NULL)
  }
  for(level in c('coarse', 'normal')){
    level_file <- file.path(sc_dir, paste0('celltype_annotations_', level, '.rds'))
    if(file.exists(level_file)){
      groups <- celltype_group_map(annotations, readRDS(level_file))
      if(!is.null(groups)){
        message('music-hc: cell types are grouped by the ', level, ' annotation level')
        return(groups)
      }
    }
  }
  stop(paste0('music-hc needs a coarser annotation level that groups the cell types ',
              '(celltype_annotations_coarse.rds or celltype_annotations_normal.rds) in ', sc_dir))
}

# Hierarchical MuSiC (music_prop.cluster): the proportions of the groups are estimated first,
# then the cell types inside each group. As in the MuSiC paper, genes with similar expression in
# the cell types of a group are used for the first step and the remaining genes for the second.
# The paper gives no threshold: for each group, the share marker_fraction of its expressed genes
# with the highest coefficient of variation between the cell-type means goes to the second step
# (parameter music_hc_marker_fraction). With very small shares the second step can fail with a
# singular matrix (seen at 0.01 and 0.02 with the lambrechts reference).
deconvolute_music_hc <- function(bulk_matrix, sc_matrix, annotations, sc_batch, celltype_groups,
                                 marker_fraction = 0.1){
  if(is.null(celltype_groups)){
    stop('music-hc: no groups of cell types were passed (argument celltype_groups)')
  }
  cell_types <- unique(as.character(annotations))
  if(!all(cell_types %in% names(celltype_groups))){
    stop(paste0('music-hc: no group for cell type(s): ', paste(setdiff(cell_types, names(celltype_groups)), collapse = ', ')))
  }
  if(is.na(marker_fraction) || marker_fraction <= 0 || marker_fraction >= 1){
    stop('music-hc: music_hc_marker_fraction has to be between 0 and 1')
  }

  # MuSiC calls Biobase functions without their namespace
  require(Biobase)

  sc_matrix <- as.matrix(sc_matrix)
  annotations <- as.character(annotations)
  pheno <- data.frame(batchId = sc_batch, cellType = annotations,
                      cellGroup = unname(celltype_groups[annotations]),
                      row.names = colnames(sc_matrix), check.names = FALSE, stringsAsFactors = FALSE)
  sc_eset <- Biobase::ExpressionSet(assayData = sc_matrix, phenoData = Biobase::AnnotatedDataFrame(pheno))
  bulk_eset <- Biobase::ExpressionSet(assayData = as.matrix(bulk_matrix))

  # MuSiC (music_prop) leaves out cell types for which it cannot compute the variance between
  # subjects, e.g. a cell type with cells from one subject only. music_prop.cluster has no such
  # check and fails on them, so the same cell types are left out here.
  basis <- MuSiC::music_basis(sc_eset, clusters = 'cellType', samples = 'batchId', select.ct = cell_types, verbose = FALSE)
  valid <- (colSums(is.na(basis$Sigma)) == 0)[cell_types] & (colSums(is.na(basis$Disgn.mtx)) == 0)[cell_types] & !is.na(basis$M.S[cell_types])
  if(!all(valid)){
    message('music-hc: left out, as MuSiC does, because the variance between subjects cannot be computed: ',
            paste(cell_types[!valid], collapse = ', '))
    cell_types <- cell_types[valid]
    sc_eset <- sc_eset[, annotations %in% cell_types]
  }
  if(length(cell_types) < 2){
    stop('music-hc: fewer than two cell types with cells from more than one subject')
  }

  # group -> cell types
  clusters_type <- split(cell_types, celltype_groups[cell_types])

  # mean relative abundance of every gene in every cell type
  theta <- basis$M.theta[, cell_types, drop = FALSE]
  bulk_genes <- rownames(bulk_matrix)[rowMeans(as.matrix(bulk_matrix)) != 0]
  theta <- theta[rownames(theta) %in% bulk_genes, , drop = FALSE]

  # genes for the cell types inside every group with more than one cell type
  select_markers <- function(fraction){
    lapply(clusters_type, function(types){
      if(length(types) < 2){
        return(NULL)
      }
      means <- theta[, types, drop = FALSE]
      variation <- apply(means, 1, sd) / rowMeans(means)
      variation <- sort(variation[is.finite(variation)], decreasing = TRUE)
      names(variation)[seq_len(ceiling(fraction * length(variation)))]
    })
  }

  # With few genes the step inside a group can fail with a singular matrix (seen with the fine
  # annotation of wu at 0.1). The share is then raised, and the share that was used is stored
  # with the result.
  fractions <- unique(c(marker_fraction, c(0.25, 0.5)[c(0.25, 0.5) > marker_fraction]))
  for(fraction in fractions){
    group_markers <- select_markers(fraction)
    for(group in names(clusters_type)){
      message('music-hc: group ', group, ' = ', paste(clusters_type[[group]], collapse = ', '),
              ' (', length(group_markers[[group]]), ' genes for the cell types inside the group)')
    }
    result <- tryCatch(
      MuSiC::music_prop.cluster(
        bulk.eset = bulk_eset,
        sc.eset = sc_eset,
        group.markers = group_markers,
        groups = 'cellGroup',
        clusters = 'cellType',
        samples = 'batchId',
        clusters.type = clusters_type,
        verbose = TRUE
      )$Est.prop.weighted.cluster,
      error = function(e) e
    )
    if(!inherits(result, 'error')){
      if(fraction != marker_fraction){
        message('music-hc: used a share of ', fraction, ' instead of ', marker_fraction)
        attr(result, 'marker_fraction') <- fraction
      }
      return(result)
    }
    message('music-hc: failed with a share of ', fraction, ': ', conditionMessage(result))
    if(!grepl('singular', conditionMessage(result))){
      stop(result)
    }
  }
  stop(result)
}

# Reference and gene weights of SECRET (Lu et al. 2023, https://github.com/anlingUA/SECRET), built
# as the paper and the vignette of the package describe:
# - marker genes of every cell type against all others with Seurat::FindAllMarkers (likelihood-ratio
#   test, genes detected in at least 25 % of the cells of one of the two groups);
# - the weight of a gene is 1 / number of cell types it is a marker of;
# - the reference is the mean expression of every cell type (SECRET::scRefer).
build_secret_reference <- function(sc_matrix, annotations, sc_batch){
  require(Biobase)
  sc_matrix <- as.matrix(sc_matrix)
  annotations <- as.character(annotations)

  seurat_obj <- Seurat::CreateSeuratObject(counts = sc_matrix)
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  Seurat::Idents(seurat_obj) <- factor(annotations)
  markers <- Seurat::FindAllMarkers(seurat_obj, test.use = 'bimod', min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  marker_frequency <- table(markers$gene)

  # scRefer collects the columns of a cell type by searching its name in the column names, which
  # mixes cell types whose names contain each other ('Monocytes', 'Monocytes classical'). It is
  # therefore called with placeholders of equal length instead of the names.
  cell_types <- unique(annotations)
  placeholders <- setNames(sprintf('celltype%04dx', seq_along(cell_types)), cell_types)
  pheno <- data.frame(sampleID = as.character(sc_batch), cellType = unname(placeholders[annotations]),
                      row.names = colnames(sc_matrix), stringsAsFactors = FALSE)
  sc_eset <- Biobase::ExpressionSet(assayData = sc_matrix, phenoData = Biobase::AnnotatedDataFrame(pheno))
  reference <- SECRET::scRefer(sc_eset, ct_var = 'cellType', sample_var = 'sampleID')
  colnames(reference) <- names(placeholders)[match(colnames(reference), placeholders)]

  genes <- intersect(rownames(reference), names(marker_frequency))
  list(reference = reference[genes, , drop = FALSE],
       weights = setNames(1 / as.numeric(marker_frequency[genes]), genes))
}

# Add the cell-type sums that the ground truth of some bulk datasets needs (e.g. all T cell
# subtypes as 'T cell'). The rules are listed in assets/celltype_aggregation.csv and applied in
# the order of the file; a cell type marked with '?' is added only if the method estimated it.
aggregate_cell_types <- function(deconvolution, bulk_name, baseDir){
  rules <- read.table(paste0(baseDir, '/assets/celltype_aggregation.csv'), sep = ',', header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  rules <- rules[sapply(strsplit(rules$bulk, ';'), function(names) bulk_name %in% names), ]
  if(nrow(rules) == 0){
    return(deconvolution)
  }
  deconvolution <- as.data.frame(deconvolution)
  for(i in seq_len(nrow(rules))){
    parts <- strsplit(rules$sum_of[i], ';')[[1]]
    optional <- grepl('\\?$', parts)
    parts <- sub('\\?$', '', parts)
    parts <- parts[!optional | parts %in% colnames(deconvolution)]
    total <- 0
    for(part in parts){
      total <- total + deconvolution[[part]]
    }
    deconvolution[[rules$cell_type[i]]] <- total
  }
  deconvolution
}

escapeCelltypesAutogenes <- function(celltype){
  celltype <- gsub("\\+", "21b2c6e87f8711ec9bf265fb9bf6ab9c", celltype)
  celltype <- gsub("-", "21b2c7567f8711ec9bf265fb9bf6ab9a", celltype)
  celltype <- gsub("\\(", "21b2c7567f8711ec9bf265fb9bf6ab9f", celltype)
  celltype <- gsub(")", "21b2c7567f8711ec9bf265fb9bf6ab9g", celltype)
  return(gsub(" ", "xxxx", celltype))
}

reEscapeCelltypesAutogenes <- function(celltype){
  celltype <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "+", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9f", "(", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9g", ")", celltype)
  return(gsub("xxxx", " ", celltype))
}

# Cell-type names of a result as they are in the reference. Methods change them: Scaden replaces
# spaces and other characters by dots, CIBERSORTx spaces by 'xxxx'. A name that matches exactly one
# cell type of the reference on its letters and digits gets that name; any other name stays.
restore_celltype_names <- function(names, reference_names){
  key <- function(x) tolower(gsub('[^A-Za-z0-9]', '', gsub('xxxx', '', x)))
  reference_names <- unique(as.character(reference_names))
  reference_keys <- key(reference_names)
  vapply(names, function(name){
    if(name %in% reference_names) return(name)
    match <- reference_names[reference_keys == key(name)]
    if(length(match) == 1) match else name
  }, character(1), USE.NAMES = FALSE)
}

compute_metrics <- function(ref_facs, deconvolution, metric = c('cor', 'rmse', 'mae', 'mape'), comparison = c('sample', 'cell_type')){
  
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

  compute_mae <- function(x, y){
    mae <- mean(abs(x - y))
    return(mae)
  }
  
  compute_mape <- function(x, y){
    mape <- mean(abs(x - y)/x)
    return(mape)
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

  } else if(metric == 'mae'){
    res_metric <- results_df %>%
      group_by(!!! syms(comparison)) %>%
      dplyr::summarize(MAP = compute_mae(estimated_frac, true_frac))

  } else if(metric == 'mape'){
    res_metric <- results_df %>%
      group_by(!!! syms(comparison)) %>%
      dplyr::summarize(MAPE = compute_mape(estimated_frac, true_frac))

  } else {
    message('Selected metric is not supported in the benchmark. Please check again.')
    stop()
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

  if(is.null(fine_annotations)){
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
  
  } else if(is.null(coarse_annotations)) {

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

  
}

# This function contains all the necessary steps to adapt omnideconv to the benchmarking workflow
signature_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, 
                                       method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path, baseDir=NULL){
  
  # path to temporary directory that is inside results directory of one unique result 
  tmp_dir_path <- paste0(res_path, '/tmp/')
  ncores <- as.numeric(ncores)

  # Signature building part
  # Dependent on which method we have 
  if (method == "cibersortx") {
    read_cibersortx_credentials(paste0(baseDir, '/cibersortx_credentials.csv'))
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(tmp_dir_path,'_input')
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input, recursive=TRUE)
    }
    cx_output <- paste0(tmp_dir_path,'_output')
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output, recursive=TRUE)
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

    library(dplyr)
    
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
      verbose = TRUE
    )$basis
    
  } else if(method == "scaden"){
    unlink(tmp_dir_path, recursive=TRUE)
    if(!dir.exists(paste0(tmp_dir_path))){
      dir.create(tmp_dir_path, recursive=TRUE)
    }
    # create output directories for scaden models, otherwise scaden might have permission issues when creating them
    #dir.create(paste0(res_path, '/model/m256'), recursive=T)
    #dir.create(paste0(res_path, '/model/m512'), recursive=T)
    #dir.create(paste0(res_path, '/model/m1024'), recursive=T)

    signature <- omnideconv::build_model_scaden(
      sc_matrix,
      annotations,
      bulk_matrix,
      temp_dir = tmp_dir_path,
      model_path = paste0(res_path, '/model'),
      verbose = TRUE
    )
  } else if (method == 'secret'){
    signature <- build_secret_reference(sc_matrix, annotations, sc_batch)

  } else if (method %in% c('autogenes', 'bayesprism', 'bisque', 'music', 'music-hc', 'rectangle')){
    signature <- NULL
      
  } else {
    message('Selected method is not supported in the benchmark. Please check again.')
    stop()
  }

  message('Finished signature building.')
  
  return(signature)
  
}

# read the CX credentials and supply them to omnideonv
read_cibersortx_credentials <- function(credentials_file){
  cred <- read.csv(credentials_file)

  omnideconv::set_cibersortx_credentials(cred$email,
                                         cred$token)  
}

# This function contains all the necessary steps to adapt omnideconv to the benchmarking workflow
deconvolution_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, 
                                           signature, method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path,
                                           rmbatch_B_mode = FALSE, rmbatch_S_mode = FALSE, baseDir=NULL,
                                           celltype_groups = NULL, method_settings = list()){
  
  # path to temporary directory that is inside results directory of one unique result 
  tmp_dir_path <- paste0(res_path, '/tmp/')
  ncores <- as.numeric(ncores)

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

    read_cibersortx_credentials(paste0(baseDir, '/cibersortx_credentials.csv')) 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(tmp_dir_path,'_input')
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input, recursive=TRUE)
    }
    cx_output <- paste0(tmp_dir_path,'_output')
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output, recursive=TRUE)
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

  } else if (method == "secret"){
    # with_unknown: also estimate the proportion of cell types that the reference does not contain.
    # The estimates of the known cell types then sum to at most 1 and the rest is 'unknown'.
    with_unknown <- is.null(method_settings$with_unknown) || toupper(method_settings$with_unknown) == 'TRUE'
    genes <- intersect(rownames(signature$reference), rownames(bulk_matrix))
    deconvolution <- SECRET::SECRET(
      bulk = bulk_matrix[genes, , drop = FALSE],
      reference = signature$reference[genes, , drop = FALSE],
      withUnknown = with_unknown,
      w = signature$weights[genes],
      yNorm = 'cpm',
      bNorm = 'cpm'
    )[[1]]

  } else if (method == "music-hc"){
    marker_fraction <- if(is.null(method_settings$marker_fraction)) 0.1 else as.numeric(method_settings$marker_fraction)
    deconvolution <- deconvolute_music_hc(bulk_matrix, sc_matrix, annotations, sc_batch, celltype_groups,
                                          marker_fraction = marker_fraction)

  } else if (method == "scdc") {    
    deconvolution <- omnideconv::deconvolute_scdc(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = annotations,
      batch_ids = sc_batch,
      verbose = TRUE,
    )

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
      temp_dir = tmp_dir_path, 
      verbose = TRUE
    )
    #unlink(tmp_dir_path, recursive=TRUE)

  } else if (method == "rectangle") {

    rp <- reticulate::import("rectanglepy")

    #ad <- anndata::AnnData(X = t(sc_matrix),
    #                       obs = data.frame('cell_type' = annotations, row.names = colnames(sc_matrix)),
    #                       var = data.frame('genes' = rownames(sc_matrix), row.names = rownames(sc_matrix)))

    deconvolution <- rp$rectangle(sc_matrix, bulks = as.data.frame(t(bulk_matrix)), n_cpus = as.integer(ncores))[[1]]

  } else {
    message('Selected method is not supported in the benchmark. Please check again.')
    stop()
  }
  
  return(deconvolution)
}

