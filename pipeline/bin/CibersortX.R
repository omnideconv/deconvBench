config_env <- new.env()

verbose_wrapper <- function(verbose) {
  return(function(method) {
    if (!verbose) {
      suppressMessages(method)
    } else {
      method
    }
  })
}

set_cibersortx_credentials <- function(email, token) {
  assign("cibersortx_email", email, envir = config_env)
  assign("cibersortx_token", token, envir = config_env)
}


build_model_cibersortx <- function(single_cell_object, cell_type_annotations,
                                   container = c("docker", "singularity"),
                                   verbose = FALSE, input_dir = NULL,
                                   output_dir = NULL, display_heatmap = FALSE,
                                   k_max = 999, ...) {
  if (is.null(single_cell_object)) {
    stop("Parameter 'single_cell_object' is missing or null, but it is required.")
  }
  if (is.null(cell_type_annotations)) {
    stop("Parameter 'cell_type_annotations' is missing or null, but it is required.")
  }
  if (container == "docker"){
    if (!docker_available()) {
      message(
        "Installation of docker can not be found. Please check whether you can ",
        "call 'docker' in the command line and get a help menu"
      )
      return(NULL)
    }
    if (!docker_connectable()) {
      message(
        "Error durching connection to docker. Please check whether you can ",
        "call 'docker ps' in the command line and get a (possibly empty) list and not an error ",
        "message"
      )
      return(NULL)
    }
  }
  check_credentials()

  temp_dir <- tempdir()

  if(container == "singularity"){
    system(paste0('singularity pull --dir ', temp_dir, ' docker://cibersortx/fractions'))
  }


  if (is.null(input_dir)) {
    input_dir <- temp_dir
  }
  if (is.null(output_dir)) {
    output_dir <- temp_dir
  }

  if (class(single_cell_object)[1] != "character") {
    transform_and_save_single_cell(single_cell_object, cell_type_annotations, input_dir, verbose)
    single_cell_object_filename <- "sample_file_for_cibersort.txt"
  } else {
    single_cell_object_filename <- single_cell_object
  }
  command_to_run <- create_docker_command(input_dir, output_dir, container,
    method = "create_sig",
    verbose = verbose,
    refsample = single_cell_object_filename, k_max = k_max
  )

  if (verbose) {
    message(command_to_run)
  }

  if (!verbose) {
    if (Sys.info()["sysname"] == "Windows") {
      message("The windows implementation requires verbose mode. It is now switched on.")
      verbose <- TRUE
    }
  }

  filebase <- paste0(
    "CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample",
    "_file_for_cibersort_inferred_refsample.bm.K", k_max
  )
  filename_sig_matrix <- paste0(filebase, ".txt")
  full_path <- paste0(input_dir, "/", filename_sig_matrix)

  if (file.exists(full_path)) {
    file.remove(full_path)
  }

  code <- system(command_to_run, ignore.stdout = !verbose, ignore.stderr = !verbose)
  if (code != 0) {
    message(paste0(
      "Something went wrong: Error code ", code, ". Please try again with ",
      "'verbose=TRUE'"
    ))
  }

  if (display_heatmap) {
    filename_heatmap <- paste0(filebase, ".pdf")
    Biobase::openPDF(normalizePath(paste0(output_dir, "/", filename_heatmap)))
  }

  sig_matrix <- verbose_wrapper(verbose)(as.data.frame(readr::read_tsv(
    paste0(output_dir, "/", filename_sig_matrix)
  )))
  rownames(sig_matrix) <- sig_matrix$NAME

  return(as.matrix.data.frame(sig_matrix[, -1]))
}

#' Deconvolute with CIBERSORTx
#'
deconvolute_cibersortx <- function(bulk_gene_expression, signature, verbose = FALSE,
                                   container = c("docker", "singularity"),
                                   input_dir = NULL, output_dir = NULL,
                                   display_extra_info = FALSE, label = "none", ...) {
  if (is.null(bulk_gene_expression)) {
    stop("Parameter 'bulk_gene_expression' is missing or null, but it is required.")
  }
  if (is.null(signature)) {
    stop("Parameter 'signature' is missing or null, but it is required.")
  }

  if (container == "docker"){
    if (!docker_available()) {
      message(
        "Installation of docker can not be found. Please check whether you can ",
        "call 'docker' in the command line and get a help menu"
      )
      return(NULL)
    }
    if (!docker_connectable()) {
      message(
        "Error durching connection to docker. Please check whether you can ",
        "call 'docker ps' in the command line and get a (possibly empty) list and not an error ",
        "message"
      )
      return(NULL)
    }
  }


  check_credentials()
  temp_dir <- tempdir()

  
  if(container == "singularity"){
    system(paste0('singularity pull --dir ', temp_dir, ' docker://cibersortx/fractions'))
  }


  if (is.null(input_dir)) {
    input_dir <- temp_dir
  }
  if (is.null(output_dir)) {
    output_dir <- temp_dir
  }

  if (class(signature)[1] != "character") {
    sig <- paste0(input_dir, "/signature_matrix.txt")
    readr::write_tsv(data.frame("NAME" = rownames(signature), signature), sig)
    sigmatrix_filename <- "signature_matrix.txt"
  } else {
    sigmatrix_filename <- signature
  }
  if (class(bulk_gene_expression)[1] != "character") {
    transform_and_save_bulk(bulk_gene_expression, input_dir, verbose)
    bulk_gene_expression_filename <- "mixture_file_for_cibersort.txt"
  } else {
    bulk_gene_expression_filename <- bulk_gene_expression
  }
  unique_id <- uuid::UUIDgenerate(TRUE)
  if (label == "none") {
    label <- unique_id
  } else {
    label <- paste0(label, "_", unique_id)
  }

  filename_cell_props <- paste0("CIBERSORTx_", label, "_Results.txt")
  cell_props_full_path <- paste0(output_dir, "/", filename_cell_props)

  command_to_run <- create_docker_command(input_dir, output_dir, container,
    method = "impute_cell_fractions", verbose = verbose,
    sigmatrix = sigmatrix_filename, mixture <- bulk_gene_expression_filename, label = label
  )
  if (verbose) {
    message(command_to_run)
  }

  if (!verbose) {
    if (Sys.info()["sysname"] == "Windows") {
      message("The windows implementation requires verbose mode. It is now switched on.")
      verbose <- TRUE
    }
  }

  code <- system(command_to_run, ignore.stdout = !verbose, ignore.stderr = !verbose)
  if (code != 0) {
    message(paste0(
      "Something went wrong: Error code ", code, ". Please try again with ",
      "'verbose=TRUE'"
    ))
  }

  cell_props_tmp <- verbose_wrapper(verbose)(as.data.frame(readr::read_tsv(
    cell_props_full_path
  )))
  cell_props <- cell_props_tmp
  rm(cell_props_tmp)
  rownames(cell_props) <- cell_props$Mixture
  cell_props <- cell_props[, -1]

  extra_cols <- c("P.value", "Correlation", "RMSE", "P-value")
  if (display_extra_info) {
    print(cell_props[, extra_cols])
  }

  cell_props <- cell_props[, !names(cell_props) %in% extra_cols]
  colnames(cell_props) <- colnames(signature)

  return(as.matrix.data.frame(cell_props))
}


transform_and_save_single_cell <- function(sc_matrix, cell_types, path, verbose = FALSE) {
  colnames(sc_matrix) <- cell_types
  output <- rbind(colnames(sc_matrix), sc_matrix)
  rownames(output) <- c("GeneSymbol", rownames(sc_matrix))
  output <- data.frame("GeneSymbol" = rownames(output), output)
  output_file <- paste0(path, "/sample_file_for_cibersort.txt")
  readr::write_tsv(output, output_file, col_names = FALSE)
  if (verbose) {
    message(paste(
      "Single cell matrix was saved successfully and can be found at: ",
      output_file
    ))
  }
  return(output_file)
}

transform_and_save_bulk <- function(bulk, path, verbose = FALSE) {
  output_file <- paste0(path, "/mixture_file_for_cibersort.txt")
  readr::write_tsv(data.frame("Gene" = rownames(bulk), bulk), output_file)
  if (verbose) {
    message(paste(
      "Bulk data matrix was saved successfully and can be found at: ",
      output_file
    ))
  }
  return(output_file)
}


create_docker_command <- function(in_dir, out_dir,
                                  container = c("docker", "singularity"),
                                  method = c("create_sig", "impute_cell_fractions"),
                                  verbose = FALSE, ...) {

  if(container=='docker'){
    base <- paste0(
      "docker run -v ", in_dir, ":/src/data:z -v ", out_dir,
      ":/src/outdir:z cibersortx/fractions --single_cell TRUE"
    )
  } else {
    base <- paste0(
      "singularity exec -c -B ",in_dir, "/:/src/data -B ", in_dir, "/:/src/outdir ", in_dir, "/fractions_latest.sif /src/CIBERSORTxFractions --single_cell TRUE"
    )
  }


  if (verbose) {
    base <- paste(base, "--verbose TRUE")
  }
  check_credentials()
  credentials <- paste(
    "--username", get("cibersortx_email", envir = config_env), "--token",
    get("cibersortx_token", envir = config_env)
  )
  return(paste(base, credentials, get_method_options(method, ...)))
}

get_method_options <- function(method = c("create_sig", "impute_cell_fractions"), ...) {
  if (method == "create_sig") {
    return(get_signature_matrix_options(...))
  } else if (method == "impute_cell_fractions") {
    return(get_cell_fractions_options(...))
  } else {
    stop(paste("Method", method, "is not valid"))
  }
}


get_signature_matrix_options <- function(refsample, g_min = 300, g_max = 500, q_value = 0.01,
                                         filter = FALSE, k_max = 999, remake = FALSE,
                                         replicates = 5, sampling = 0.5, fraction = 0.75) {
  return(paste(
    "--refsample", refsample, "--G.min", g_min, "--G.max", g_max, "--q.value", q_value, "--filter",
    filter, "--k.max", k_max, "--remake", remake, "--replicates", replicates, "--sampling",
    sampling, "--fraction", fraction
  ))
}



get_cell_fractions_options <- function(sigmatrix, mixture, perm = 0, label = "none",
                                       rmbatch_B_mode = FALSE, rmbatch_S_mode = FALSE,
                                       source_GEPs = sigmatrix, qn = FALSE,
                                       absolute = FALSE, abs_method = "sig.score") {
  return(paste(
    "--mixture", mixture, "--sigmatrix", sigmatrix, "--perm", perm, "--label", label,
    "--rmbatchBmode", rmbatch_B_mode, "--rmbatchSmode", rmbatch_S_mode, "--sourceGEPs", source_GEPs,
    "QN", qn, "--absolute", absolute, "--abs_method", abs_method
  ))
}

#' Checks that the email and token variables are set
#'
check_credentials <- function() {
  assertthat::assert_that(exists("cibersortx_email", envir = config_env),
    msg = paste(
      "CIBERSORTx email for credentials is missing. Please call",
      "set_cibersortx_credentials(email,token) first."
    )
  )
  assertthat::assert_that(exists("cibersortx_token", envir = config_env),
    msg = paste(
      "CIBERSORTx token for credentials is missing. Please call",
      "set_cibersortx_credentials(email,token) first."
    )
  )
}
