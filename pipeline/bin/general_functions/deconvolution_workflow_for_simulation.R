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



# This function contains all the necessary steps to adapt omnideconv to the benchmarking cluster
signature_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path, spillover = FALSE){
  
  # Signature building part
  # Dependent on which method we have 
  if(method == 'autogenes' | method == 'scaden'){
    annotations <- escapeCelltypesAutogenes(annotations)
    
    signature_dir <- paste0(res_path, '/tmp/', method, '_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot/")
    
    unlink(signature_dir, recursive=TRUE)
    if(!dir.exists(paste0(signature_dir))){
      dir.create(signature_dir)
    }
    
    if(method == 'autogenes'){
      signature <- omnideconv::build_model_autogenes(
        sc_matrix,
        annotations,
        output_dir = signature_dir,
        verbose = TRUE
      )
    } else {
      signature <- omnideconv::build_model_scaden(
        sc_matrix,
        annotations,
        bulk_matrix,
        temp_dir = scaden_tmp,
        verbose = TRUE
      )
    }
    
    old_signatures <- list.files(res_path, ".pickle", full.names = TRUE) 
    
    if(length(old_signatures > 0)){
      print('Found old signatures in output directory; they will be removed.')
      sapply(old_signatures, file.remove)
    }
    
    #We copy the new signature to the results directory
    file.copy(signature, res_path, recursive = TRUE)
    
    signature <- list.files(res_path, ".pickle", full.names = TRUE)
    
  } else if(method == 'cibersortx'){
    
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                           " 721a387e91c495174066462484674cb8") 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(res_path, '/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    
    cx_output <- paste0(res_path, '/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output)
    }
    
    signature <- omnideconv::build_model_cibersortx(
      sc_matrix, 
      sc_celltype_annotations, 
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
      ncores = ncores
    )
    
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
  
  saveRDS(signature, file=paste0(res_path, "/signature.rds"))
  
  signature
  
}


deconvolution_workflow_general <- function(sc_matrix, annotations, annotation_category, sc_ds, sc_norm, sc_batch, signature, method, bulk_matrix, bulk_name, bulk_norm, ncores, res_path){
  
  if(method == 'autogenes'){
    
    signature_dir <- paste0(res_path, '/tmp/autogenes_tmp_',sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot/")
    
    deconvolution <- omnideconv::deconvolute_autogenes(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      max_iter = 1000000
    )$proportions
    colnames(deconvolution) <- reEscapeCelltypesAutogenes(colnames(deconvolution))
    unlink(signature_dir, recursive = TRUE)
    file.remove(signature)
    
  } else if(method == 'scaden'){
    
    scaden_tmp <- paste0(res_path, '/tmp/scaden_tmp_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    deconvolution <- omnideconv::deconvolute_scaden(
      bulk_gene_expression = bulk_matrix, 
      signature = signature,
      temp_dir = scaden_tmp
    )
    unlink(scaden_tmp, recursive=TRUE)
    
  } else if(method == 'cibersortx'){
    
    omnideconv::set_cibersortx_credentials("lorenzo.merotto@studenti.unipd.it",
                                           " 721a387e91c495174066462484674cb8") 
    # CibersortX does not work with tmp directories in a Docker in Docker setup
    # --> created fixed input and output directories!
    cx_input <- paste0(res_path, '/tmp/cibersortx_input_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_input))){
      dir.create(cx_input)
    }
    
    cx_output <- paste0(res_path, '/tmp/cibersortx_output_', sc_ds, "_", sc_norm, "_", bulk_name, "_", bulk_norm, "_", annotation_category, "_annot")
    
    if(!dir.exists(paste0(cx_output))){
      dir.create(cx_output)
    }
    
    deconvolution <- omnideconv::deconvolute_cibersortx(
      bulk_gene_expression = bulk_matrix, 
      signature = signature, 
      container = 'docker',
      verbose = TRUE,
      input_dir = cx_input,
      output_dir = cx_output
    )
    unlink(cx_input, recursive=TRUE)
    
  } else if(method == 'dwls'){
    
    deconvolution <- NULL
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
        colnames(deconvolution) <<- as.character(colnames(signature))
        rownames(deconvolution) <<- as.character(colnames(bulk_matrix))
      }
    )
    
  } else if(method == 'bayesprism'){
    
    deconvolution <- omnideconv::deconvolute_bayesprism(
      bulk_gene_expression = bulk_matrix, 
      single_cell_object = sc_matrix, 
      cell_type_annotations = annotations, 
      n_cores = ncores,
      species = 'hs'
    )$theta
    
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
  
  colnames(deconvolution) <- gsub('\\.',' ',colnames(deconvolution))
  deconvolution
}

