# Checks of a CI run (tests/ci/run_ci.sh). Fails with a list of problems; warnings are printed only.
#
# Usage: Rscript check_results.R <suite> <run_dir> <pipeline_dir> [--update-expected]
#   --update-expected   write the MuSiC results of this run as the new expected values

args <- commandArgs(trailingOnly = TRUE)
suite <- args[1]; run_dir <- args[2]; pipeline_dir <- args[3]
update_expected <- '--update-expected' %in% args
ci_dir <- file.path(pipeline_dir, 'tests', 'ci')
data_dir <- file.path(dirname(pipeline_dir), 'data', 'ci')

problems <- character(); warnings <- character()
fail <- function(...) problems <<- c(problems, paste0(...))
warn <- function(...) warnings <<- c(warnings, paste0(...))
check <- function(ok, ...) if (!isTRUE(ok)) fail(...)

# cell types of the CI references at every level; methods change names (dots, underscores, 'xxxx'),
# so names are compared on letters and digits
key <- function(x) tolower(gsub('[^A-Za-z0-9]', '', gsub('xxxx', '', x)))
ref_types <- unique(unlist(lapply(list.files(file.path(data_dir, 'single_cell'), full.names = TRUE), function(d)
  unlist(lapply(list.files(d, '^celltype_annotations.*rds$', full.names = TRUE), function(f) as.character(readRDS(f)))))))
allowed <- c(ref_types, 'T cell', 'unknown', 'Unknown')

estimate <- function(x) {
  if (is.list(x) && !is.data.frame(x)) x <- x[[intersect(c('deconvolution', 'deconv.results'), names(x))[1]]]
  as.matrix(x)
}

# every result file: complete estimates whose cell types are those of the reference
check_estimates <- function(files, id) {
  check(length(files) > 0, id, ': no result files')
  for (f in files) {
    rel <- sub(paste0(run_dir, '/'), '', f, fixed = TRUE)
    m <- tryCatch(estimate(readRDS(f)), error = function(e) NULL)
    if (is.null(m) || length(m) == 0) { fail(rel, ': no estimates'); next }
    check(is.numeric(m) && !anyNA(m), rel, ': missing or non-numeric estimates')
    # sums for ground truths such as 'T cell' of hoek (see assets/celltype_aggregation.csv) do not count
    m_sum <- m[, setdiff(colnames(m), 'T cell'), drop = FALSE]
    check(all(abs(rowSums(m_sum) - 1) < 1e-3), rel, ': fractions do not sum to 1 (', paste(round(range(rowSums(m_sum)), 3), collapse = '-'), ')')
    unknown_names <- colnames(m)[!key(colnames(m)) %in% key(allowed)]
    check(length(unknown_names) == 0, rel, ': cell types not in the reference: ', paste(unknown_names, collapse = ', '))
    changed <- colnames(m)[key(colnames(m)) %in% key(allowed) & !colnames(m) %in% allowed]
    if (length(changed) > 0) warn(rel, ': cell-type names changed by the method or the pipeline: ', paste(unique(changed), collapse = ', '))
  }
}

results <- function(id, pattern = '^deconvolution.*\\.rds$') list.files(file.path(run_dir, id, 'results'), pattern, recursive = TRUE, full.names = TRUE)

# expected values: MuSiC is deterministic, so its results are compared with stored ones
compare_expected <- function(file, name) {
  m <- estimate(readRDS(file))
  m <- m[order(rownames(m)), order(colnames(m)), drop = FALSE]
  expected_file <- file.path(ci_dir, 'expected', paste0(name, '.csv'))
  if (update_expected) {
    dir.create(dirname(expected_file), showWarnings = FALSE)
    write.csv(round(m, 10), expected_file)
    return(invisible())
  }
  if (!file.exists(expected_file)) { fail('no expected values ', basename(expected_file), ' (run with --update-expected)'); return(invisible()) }
  e <- as.matrix(read.csv(expected_file, row.names = 1, check.names = FALSE))
  same <- identical(dim(m), dim(e)) && identical(rownames(m), rownames(e)) && identical(colnames(m), colnames(e)) && max(abs(m - e)) < 1e-6
  check(same, name, ': MuSiC results differ from tests/ci/expected/', basename(expected_file),
        ' (if the change is intended: run_ci.sh workflows, then check_results.R with --update-expected)')
}

if (suite == 'workflows') {
  ids <- list.dirs(run_dir, recursive = FALSE, full.names = FALSE)
  ids <- setdiff(ids, c('work', 'unit'))
  for (id in ids) check_estimates(results(id), id)

  # main: metrics only for the bulk with ground truth
  general <- file.path(run_dir, 'main', 'results', 'results_general')
  check(file.exists(file.path(general, 'music_sc-ci_counts_bulk-ci_tpm_ct0_rep0', 'results_metric.rds')), 'main: no metrics for bulk-ci')
  check(file.exists(file.path(general, 'music_sc-ci_counts_bulk-nofacs_tpm_ct0_rep0', 'deconvolution.rds')), 'main: bulk-nofacs was not deconvolved')
  check(!file.exists(file.path(general, 'music_sc-ci_counts_bulk-nofacs_tpm_ct0_rep0', 'results_metric.rds')), 'main: metrics for bulk-nofacs, which has no ground truth')
  compare_expected(file.path(general, 'music_sc-ci_counts_bulk-ci_tpm_ct0_rep0', 'deconvolution.rds'), 'main_music_bulk-ci')

  # unknown content: the pseudo-bulks contain only the known cell types and the unknown one
  known <- c('B cells', 'NK cells', 'T cells CD4 conv', 'Monocytes'); unknown_type <- 'Platelets.CD41+'
  for (id in grep('^simulation_unknown_content', ids, value = TRUE)) {
    for (f in list.files(file.path(run_dir, id, 'input', 'preprocess'), '_facs\\.rds$', recursive = TRUE, full.names = TRUE)) {
      facs <- as.matrix(readRDS(f))
      extra <- setdiff(rownames(facs)[rowSums(facs) > 0], c(known, unknown_type))
      check(length(extra) == 0, id, ': pseudo-bulk ', basename(f), ' contains other cell types: ', paste(extra, collapse = ', '))
      fraction <- as.numeric(sub('.*_([0-9.]+)_facs\\.rds$', '\\1', f))
      check(all(abs(facs[unknown_type, ] - fraction) < 0.05), id, ': ', basename(f), ': unknown fraction is not ', fraction)
    }
  }
  compare_expected(results('simulation_unknown_content', '_0.5\\.rds$')[1], 'unknown_content_music_0.5')

  # missing cell types: every combination, and the hoek ground truth sums the T cell subtypes
  for (id in grep('^impact_missing_cell_types', ids, value = TRUE)) {
    res <- results(id)
    check(length(res) == 6, id, ': ', length(res), ' results instead of 6 (2 bulks x 3 references)')
    for (f in res[grepl('_hoek_', basename(dirname(res)))]) check('T cell' %in% colnames(estimate(readRDS(f))), id, ': no T cell sum in ', basename(dirname(f)))
  }

  # resolution: three levels with coarse annotation; ground truth of every level in the annotation
  normal <- as.character(readRDS(file.path(data_dir, 'single_cell', 'sc-ci', 'celltype_annotations_normal.rds')))
  coarse <- as.character(readRDS(file.path(data_dir, 'single_cell', 'sc-ci', 'celltype_annotations_coarse.rds')))
  for (id in grep('^simulation_resolution_analysis', ids, value = TRUE)) {
    res <- results(id)
    for (level in c('fine', 'normal', 'coarse')) check(any(grepl(paste0('_', level, '_annot/'), res)), id, ': no results at the ', level, ' level')
    for (f in grep('_(normal|coarse)_annot/', res, value = TRUE)) {
      truth <- rownames(readRDS(f)$true_cell_fractions)
      level_types <- if (grepl('_normal_annot/', f)) normal else coarse
      check(all(truth %in% level_types), id, ': ground truth of ', sub('.*(normal|coarse)_annot.*', '\\1', f),
            ' level has cell types that are not in the annotation (cell_type_mappings.csv?): ', paste(setdiff(truth, level_types), collapse = ', '))
    }
  }
  compare_expected(grep('_fine_annot/', results('simulation_resolution_analysis'), value = TRUE)[1], 'resolution_music_fine')

  # reference splits
  splits <- list.dirs(file.path(run_dir, ids, 'input', 'preprocess', 'reference_split'), recursive = FALSE)
  check(length(splits) > 0, 'no reference splits were written')
  for (s in splits) {
    a <- read.csv(file.path(s, 'assignment.csv'))
    parts <- tapply(a$part, a$donor, function(p) length(unique(p[p != 'not used'])))
    if (grepl('_bycell$', s)) {
      fine <- as.character(readRDS(file.path(data_dir, 'single_cell', 'sc-ci', 'celltype_annotations_fine.rds')))
      both <- tapply(a$part, fine[a$cell], function(p) all(c('simulation', 'deconvolution') %in% p))
      check(all(both), basename(s), ': cell types not in both parts: ', paste(names(both)[!both], collapse = ', '))
      check(any(parts > 1), basename(s), ': no donor is in both parts')
    } else {
      check(all(parts == 1), basename(s), ': donors in both parts: ', paste(names(parts)[parts > 1], collapse = ', '))
    }
  }
  # the same parameters give the same split in every workflow
  for (kind in c('_bycell$', '_cells1$')) {
    files <- file.path(grep(kind, splits, value = TRUE), 'assignment.csv')
    if (length(files) > 1) check(length(unique(tools::md5sum(files))) == 1, 'splits ', kind, ' differ between workflows')
  }
} else {
  methods <- switch(suite, standard = c('bayesprism', 'bisque', 'dwls', 'music', 'music-hc', 'scaden', 'scdc'), suite)
  general <- file.path(run_dir, paste0('main-', suite), 'results', 'results_general')
  check_estimates(results(paste0('main-', suite)), paste0('main-', suite))
  # MuSiC and music-hc leave out the cell type whose cells come from a single donor
  single_donor <- 'pDC_c4-IRF7+'
  for (m in methods) {
    for (bulk in c('bulk-ci', 'bulk-nofacs')) {
      dir <- list.files(general, paste0('^', m, '_sc-ci_[a-z]+_', bulk, '_'), full.names = TRUE)
      if (length(dir) != 1) { fail(m, ' / ', bulk, ': no result folder'); next }
      est <- estimate(readRDS(file.path(dir, 'deconvolution.rds')))
      expected_types <- setdiff(normal_types <- as.character(readRDS(file.path(data_dir, 'single_cell', 'sc-ci', 'celltype_annotations_normal.rds'))), if (m %in% c('music', 'music-hc')) single_donor)
      missing_types <- unique(expected_types[!key(expected_types) %in% key(colnames(est))])
      check(length(missing_types) == 0, m, ' / ', bulk, ': cell types missing in the result: ', paste(missing_types, collapse = ', '))
      check(file.exists(file.path(dir, 'results_metric.rds')) == (bulk == 'bulk-ci'), m, ' / ', bulk, ': metrics ', if (bulk == 'bulk-ci') 'missing' else 'computed without ground truth')
    }
  }
}

if (length(warnings) > 0) { cat('Warnings:\n'); cat(paste0('  ', unique(warnings), '\n'), sep = '') }
if (update_expected) cat('Expected values written to', file.path(ci_dir, 'expected'), '\n')
if (length(problems) > 0) {
  cat('Problems:\n'); cat(paste0('  ', problems, '\n'), sep = '')
  quit(status = 1)
}
cat('All checks passed (', suite, ')\n', sep = '')
