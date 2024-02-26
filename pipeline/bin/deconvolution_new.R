
library(omnideconv)

escape_special_chars <- function(string) {
  if (is.null(string)) {
    return(NULL)
  }
  string <- gsub("\u0020", "21b29fb07f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # " "
  string <- gsub("\u0021", "21b29fb17f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "!"
  string <- gsub("\u0022", "21b29fb27f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # """
  string <- gsub("\u0023", "21b29fb37f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "#"
  string <- gsub("\u0024", "21b29fb47f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "$"
  string <- gsub("\u0025", "21b29fb57f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "%"
  string <- gsub("\u0026", "21b29fb67f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "&"
  string <- gsub("\u0027", "21b29fb77f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "'"
  string <- gsub("\u0028", "21b29fb87f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "("
  string <- gsub("\u0029", "21b29fb97f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ")"
  string <- gsub("\u002A", "21b29fba7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "*"
  string <- gsub("\u002B", "21b2c6e87f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "+"
  string <- gsub("\u002C", "21b2c6e97f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ","
  string <- gsub("\u002D", "21b2c7567f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "-"
  string <- gsub("\u002E", "21b2c7577f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "."
  string <- gsub("\u002F", "21b2c7587f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "/"
  string <- gsub("\u003A", "21b2c7597f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ":"
  string <- gsub("\u003B", "21b2c75a7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ";"
  string <- gsub("\u003C", "21b2c75b7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "<"
  string <- gsub("\u003D", "21b2c75c7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "="
  string <- gsub("\u003E", "21b2c75d7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # ">"
  string <- gsub("\u003F", "21b2c75e7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "?"
  string <- gsub("\u0040", "21b2c75f7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "@"
  string <- gsub("\u005B", "21b2c7607f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "["
  string <- gsub("\u005C", "21b2ee347f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "\"
  string <- gsub("\u005D", "21b2ee357f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "]"
  string <- gsub("\u005E", "21b2ee367f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "^"
  string <- gsub("\u005F", "21b2ee377f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "_"
  string <- gsub("\u0060", "21b2ee387f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "`"
  string <- gsub("\u00B4", "21b2ee397f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "´"
  string <- gsub("\u007B", "21b2ee3a7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "{"
  string <- gsub("\u007C", "21b2ee3b7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "|"
  string <- gsub("\u007D", "21b2ee3c7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "}"
  string <- gsub("\u007E", "21b2ee3d7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "~"
  string <- gsub("\u00A7", "21b2ee3e7f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "§"
  string <- gsub("\u00DF", "21b315447f8711ec9bf265fb9bf6ab9c", string, fixed = TRUE) # "ß"

  return(string)
}


#' Removes the substitutions and turns them back into the special characters
#'
#' @param string The string to be de-escaped
#'
#' @return The String with special characters
#'
deescape_special_chars <- function(string) {
  if (is.null(string)) {
    return(NULL)
  }
  string <- gsub("21b29fb07f8711ec9bf265fb9bf6ab9c", "\u0020", string, fixed = TRUE) # " "
  string <- gsub("21b29fb17f8711ec9bf265fb9bf6ab9c", "\u0021", string, fixed = TRUE) # "!"
  string <- gsub("21b29fb27f8711ec9bf265fb9bf6ab9c", "\u0022", string, fixed = TRUE) # """
  string <- gsub("21b29fb37f8711ec9bf265fb9bf6ab9c", "\u0023", string, fixed = TRUE) # "#"
  string <- gsub("21b29fb47f8711ec9bf265fb9bf6ab9c", "\u0024", string, fixed = TRUE) # "$"
  string <- gsub("21b29fb57f8711ec9bf265fb9bf6ab9c", "\u0025", string, fixed = TRUE) # "%"
  string <- gsub("21b29fb67f8711ec9bf265fb9bf6ab9c", "\u0026", string, fixed = TRUE) # "&"
  string <- gsub("21b29fb77f8711ec9bf265fb9bf6ab9c", "\u0027", string, fixed = TRUE) # "'"
  string <- gsub("21b29fb87f8711ec9bf265fb9bf6ab9c", "\u0028", string, fixed = TRUE) # "("
  string <- gsub("21b29fb97f8711ec9bf265fb9bf6ab9c", "\u0029", string, fixed = TRUE) # ")"
  string <- gsub("21b29fba7f8711ec9bf265fb9bf6ab9c", "\u002A", string, fixed = TRUE) # "*"
  string <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "\u002B", string, fixed = TRUE) # "+"
  string <- gsub("21b2c6e97f8711ec9bf265fb9bf6ab9c", "\u002C", string, fixed = TRUE) # ","
  string <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9c", "\u002D", string, fixed = TRUE) # "-"
  string <- gsub("21b2c7577f8711ec9bf265fb9bf6ab9c", "\u002E", string, fixed = TRUE) # "."
  string <- gsub("21b2c7587f8711ec9bf265fb9bf6ab9c", "\u002F", string, fixed = TRUE) # "/"
  string <- gsub("21b2c7597f8711ec9bf265fb9bf6ab9c", "\u003A", string, fixed = TRUE) # ":"
  string <- gsub("21b2c75a7f8711ec9bf265fb9bf6ab9c", "\u003B", string, fixed = TRUE) # ";"
  string <- gsub("21b2c75b7f8711ec9bf265fb9bf6ab9c", "\u003C", string, fixed = TRUE) # "<"
  string <- gsub("21b2c75c7f8711ec9bf265fb9bf6ab9c", "\u003D", string, fixed = TRUE) # "="
  string <- gsub("21b2c75d7f8711ec9bf265fb9bf6ab9c", "\u003E", string, fixed = TRUE) # ">"
  string <- gsub("21b2c75e7f8711ec9bf265fb9bf6ab9c", "\u003F", string, fixed = TRUE) # "?"
  string <- gsub("21b2c75f7f8711ec9bf265fb9bf6ab9c", "\u0040", string, fixed = TRUE) # "@"
  string <- gsub("21b2c7607f8711ec9bf265fb9bf6ab9c", "\u005B", string, fixed = TRUE) # "["
  string <- gsub("21b2ee347f8711ec9bf265fb9bf6ab9c", "\u005C", string, fixed = TRUE) # "\"
  string <- gsub("21b2ee357f8711ec9bf265fb9bf6ab9c", "\u005D", string, fixed = TRUE) # "]"
  string <- gsub("21b2ee367f8711ec9bf265fb9bf6ab9c", "\u005E", string, fixed = TRUE) # "^"
  string <- gsub("21b2ee377f8711ec9bf265fb9bf6ab9c", "\u005F", string, fixed = TRUE) # "_"
  string <- gsub("21b2ee387f8711ec9bf265fb9bf6ab9c", "\u0060", string, fixed = TRUE) # "`"
  string <- gsub("21b2ee397f8711ec9bf265fb9bf6ab9c", "\u00B4", string, fixed = TRUE) # "´"
  string <- gsub("21b2ee3a7f8711ec9bf265fb9bf6ab9c", "\u007B", string, fixed = TRUE) # "{"
  string <- gsub("21b2ee3b7f8711ec9bf265fb9bf6ab9c", "\u007C", string, fixed = TRUE) # "|"
  string <- gsub("21b2ee3c7f8711ec9bf265fb9bf6ab9c", "\u007D", string, fixed = TRUE) # "}"
  string <- gsub("21b2ee3d7f8711ec9bf265fb9bf6ab9c", "\u007E", string, fixed = TRUE) # "~"
  string <- gsub("21b2ee3e7f8711ec9bf265fb9bf6ab9c", "\u00A7", string, fixed = TRUE) # "§"
  string <- gsub("21b315447f8711ec9bf265fb9bf6ab9c", "\u00DF", string, fixed = TRUE) # "ß"
  return(string)
}


convert_to_matrix <- function(object, cell_type_annotations, cell_type_column_name = NULL) {
  if (!is.null(object)) {
    if (class(object)[[1]] == "AnnDataR6") {
      object <- anndata_to_singlecellexperiment(object)
    }

    if (class(object)[[1]] == "SingleCellExperiment") {
      matrix_and_annotation <-
        singlecellexperiment_to_matrix(object,
          cell_type_column_name = cell_type_column_name
        )
      object <- matrix_and_annotation$matrix
      if (is.null(cell_type_annotations)) {
        if (is.null(cell_type_column_name)) {
          stop(
            "Either provide cell type annotations as vector (cell_type_annotations) or the ",
            "name of the column that stores label information!"
          )
        } else {
          cell_type_annotations <- matrix_and_annotation$annotation_vector
        }
      }
    }

    if (class(object)[[1]] != "matrix") {
      object <- as.matrix(object)
    }
  }

  return(list(matrix = object, cell_type_annotations = cell_type_annotations))
}


#' Deconvolution
#'
#' @param bulk_gene_expression A matrix or dataframe with the bulk data. Rows are genes, columns
#'   are samples.
#' @param signature The signature matrix.
#' @param method A string specifying the method.
#'   Supported methods are 'bisque', 'momf', 'dwls', 'scaden', 'cibersortx' and 'autogenes'
#' @param single_cell_object Needed for deconvolution with MOMF and Bisque. Defaults to NULL.
#'   Alternatively a SingleCellExperiment or an AnnData object can be provided. In that case, note
#'   that cell-type labels need to be indicated either directly providing a vector
#'   (cell_type_annotations) or by indicating the column name that indicates the cell-type labels
#'   (cell_type_column_name). (Anndata: obs object, SingleCellExperiment: colData object)
#' @param cell_type_annotations Needed for deconvolution with Bisque, MuSiC and SCDC.
#'   Defaults to NULL.
#' @param batch_ids A vector of the ids of the samples or individuals. Defaults to NULL.
#' @param verbose Whether to produce an output on the console.
#' @param ... Additional parameters, passed to the algorithm used.
#' @param cell_type_column_name Name of the column in (Anndata: obs, SingleCellExperiment: colData),
#'   that contains the cell-type labels. Is only used if no cell_type_annotations vector
#'   is provided.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are
#' individuals, columns are cell types.
#' @export
#'
#' @examples
#' # More examples can be found in the unit tests at tests/testthat/test-c-deconvolute.R
#' data("single_cell_data_1")
#' data("cell_type_annotations_1")
#' data("batch_ids_1")
#' data("bulk")
#'
#' single_cell_data <- single_cell_data_1[1:2000, 1:500]
#' cell_type_annotations <- cell_type_annotations_1[1:500]
#' batch_ids <- batch_ids_1[1:500]
#' bulk <- bulk[1:2000, ]
#'
#' signature_matrix_autogenes <- build_model(
#'   single_cell_data, cell_type_annotations, "autogenes",
#'   batch_ids
#' )
#' deconv_autogenes <- deconvolute(
#'   bulk, signature_matrix_autogenes, "autogenes"
#' )
#'
#' deconv_bisque <- deconvolute(
#'   bulk, NULL, "bisque", single_cell_data,
#'   cell_type_annotations, batch_ids
#' )
deconvolute_new <- function(bulk_gene_expression, signature, method = deconvolution_methods,
                        single_cell_object = NULL, cell_type_annotations = NULL, batch_ids = NULL,
                        cell_type_column_name = NULL, verbose = FALSE, ...) {
  if (length(method) > 1) {
    stop(
      "Please only specify one method and not ", length(method), ": ",
      paste(method, collapse = ", ")
    )
  }
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  method <- tolower(method)
  

  # Converting all other data types into a matrix
  matrix_and_annotation <- convert_to_matrix(
    single_cell_object, cell_type_annotations,
    cell_type_column_name
  )
  single_cell_object <- matrix_and_annotation$matrix
  cell_type_annotations <- matrix_and_annotation$cell_type_annotations


  # Converting all other data types into a matrix
  bulk_gene_expression <- convert_to_matrix(bulk_gene_expression, "bulk")$matrix

  # Check the input data for problems like different numbers of cells in the object and the
  # annotation or strings in the data
  #check_data(single_cell_object, cell_type_annotations, bulk_gene_expression)

  rownames(bulk_gene_expression) <- escape_special_chars(rownames(bulk_gene_expression))
  colnames(bulk_gene_expression) <- escape_special_chars(colnames(bulk_gene_expression))
  # Only do if it is a matrix or dataframe
  if ("matrix" %in% class(signature) || "data.frame" %in% class(signature)) {
    rownames(signature) <- escape_special_chars(rownames(signature))
    colnames(signature) <- escape_special_chars(colnames(signature))
  }
  if (!is.null(single_cell_object)) {
    if ("matrix" %in% class(single_cell_object) || "data.frame" %in% class(single_cell_object)) {
      rownames(single_cell_object) <- escape_special_chars(rownames(single_cell_object))
      colnames(single_cell_object) <- escape_special_chars(colnames(single_cell_object))
    } else if ("list" %in% class(single_cell_object)) {
      single_cell_object <- lapply(single_cell_object, function(sc) {
        rownames(sc) <- escape_special_chars(rownames(sc))
        colnames(sc) <- escape_special_chars(colnames(sc))
        sc
      })
    }
  }
  if (!is.null(cell_type_annotations)) {
    if ("character" %in% class(cell_type_annotations)) {
      cell_type_annotations <- escape_special_chars(cell_type_annotations)
    } else if ("list" %in% class(cell_type_annotations)) {
      cell_type_annotations <- lapply(cell_type_annotations, escape_special_chars)
    }
  }

  if (!is.null(batch_ids)) {
    if ("character" %in% class(batch_ids)) {
      batch_ids <- escape_special_chars(batch_ids)
    } else if ("list" %in% class(batch_ids)) {
      batch_ids <- lapply(batch_ids, escape_special_chars)
    }
  }

  if (verbose && method %in% c("bisque", "music", "scdc", "cpm", "cdseq") && !is.null(signature)) {
    message(
      "A signature was provided, even though you chose a method that does not use ",
      "an external one."
    )
  }

  deconv <- switch(method,
    bisque = t(deconvolute_bisque(bulk_gene_expression, single_cell_object, cell_type_annotations,
      batch_ids,
      verbose = verbose, ...
    )$bulk.props),
    momf = deconvolute_momf(bulk_gene_expression, signature, single_cell_object,
      verbose = verbose, ...
    )$cell.prop,
    scaden = deconvolute_scaden(signature, bulk_gene_expression, verbose = verbose, ...),
    dwls = deconvolute_dwls(bulk_gene_expression, signature, verbose = verbose, ...),
    cibersortx = deconvolute_cibersortx(bulk_gene_expression, signature, verbose = verbose, ...),
    autogenes = deconvolute_autogenes(bulk_gene_expression, signature,
      verbose = verbose, ...
    )$proportions,
    music = deconvolute_music(bulk_gene_expression, single_cell_object, cell_type_annotations,
      batch_ids,
      verbose = verbose, ...
    )$Est.prop.weighted,
    scdc = {
      res <- deconvolute_scdc(bulk_gene_expression, single_cell_object, cell_type_annotations,
        batch_ids,
        verbose = verbose, ...
      )
      if ("prop.est.mvw" %in% names(res)) {
        res$prop.est.mvw
      } else if ("w_table" %in% names(res)) {
        SCDC::wt_prop(res$w_table, res$prop.only)
      } else {
        message(
          "There seems to be an error, as the result of deconvolute_scdc did not ",
          "contain prop.est.mvw or w_table"
        )
        res
      }
    },
    cpm = deconvolute_cpm(bulk_gene_expression, single_cell_object, cell_type_annotations,
      verbose = verbose, ...
    )$cellTypePredictions,
    bseqsc = t(deconvolute_bseqsc(bulk_gene_expression, signature,
      verbose = verbose, ...
    )$coefficients),
    cdseq = t(deconvolute_cdseq(bulk_gene_expression, single_cell_object, cell_type_annotations,
      batch_ids,
      verbose = verbose, ...
    )$cdseq_prop_merged)
  )

  if (!is.null(deconv)) {
    # Normalize the results
    #deconv <- normalize_deconv_results(deconv)
    # Alphabetical order of celltypes
    rownames(deconv) <- deescape_special_chars(rownames(deconv))
    colnames(deconv) <- deescape_special_chars(colnames(deconv))
    deconv <- deconv[, order(colnames(deconv)), drop = FALSE]
  }
  return(deconv)
}
