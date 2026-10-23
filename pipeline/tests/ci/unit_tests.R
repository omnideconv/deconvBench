# Unit tests of pipeline helpers with the CI data. Sections that need a package of another image
# (SECRET) run only where it is installed.
#
# Usage: Rscript unit_tests.R <repo_dir> <tmp_dir>

args <- commandArgs(trailingOnly = TRUE)
repo_dir <- args[1]; tmp_dir <- args[2]
pipeline_dir <- file.path(repo_dir, 'pipeline')
sc_dir <- file.path(repo_dir, 'data', 'ci', 'single_cell', 'sc-ci')
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
suppressMessages(source(file.path(pipeline_dir, 'bin', 'utils.R')))
source(file.path(pipeline_dir, 'tests', 'ci', 'check_mappings.R'))

problems <- character()
test <- function(name, expr) {
  result <- tryCatch(expr, error = function(e) paste('error:', conditionMessage(e)))
  if (isTRUE(result)) cat('ok   ', name, '\n') else { cat('FAIL ', name, if (is.character(result)) paste0(': ', result), '\n'); problems <<- c(problems, name) }
}

counts <- readRDS(file.path(sc_dir, 'matrix_counts.rds'))
annot <- readRDS(file.path(sc_dir, 'celltype_annotations.rds'))
batch <- readRDS(file.path(sc_dir, 'batch.rds'))
coarse <- readRDS(file.path(sc_dir, 'celltype_annotations_coarse.rds'))
bulk <- as.matrix(readRDS(file.path(repo_dir, 'data', 'ci', 'bulks', 'bulk-ci', 'bulk-ci_tpm.rds')))

test('the CI data keep their edge cases', {
  is.factor(annot) && 'unused type' %in% levels(annot) && !'unused type' %in% annot &&
    length(unique(batch[annot == 'pDC_c4-IRF7+'])) == 1 &&
    all(c('Monocytes', 'Monocytes classical', 'Platelets.CD41+') %in% annot)
})

test('cell_type_mappings.csv of the CI data matches the annotations', {
  length(unlist(lapply(file.path(repo_dir, 'data', 'ci', 'single_cell', c('sc-ci', 'sc-ci-nc')), check_mapping))) == 0
})

test('check_mapping finds a misspelled fine cell type', {
  bad <- file.path(tmp_dir, 'bad_mapping'); dir.create(bad, showWarnings = FALSE)
  for (f in list.files(sc_dir, '^celltype_annotations_(fine|normal|coarse)\\.rds$')) file.copy(file.path(sc_dir, f), bad, overwrite = TRUE)
  m <- read.csv(file.path(sc_dir, 'cell_type_mappings.csv'), check.names = FALSE)
  m$Fine[m$Fine == 'Monocytes classical'] <- 'Monocyte classical'
  write.csv(m, file.path(bad, 'cell_type_mappings.csv'), row.names = FALSE)
  length(check_mapping(bad)) > 0
})

test('celltype_group_map groups by the coarse level and rejects a level that does not group', {
  groups <- celltype_group_map(as.character(annot), coarse)
  !is.null(groups) && groups[['Monocytes classical']] == groups[['Monocytes']] &&
    is.null(celltype_group_map(as.character(annot), as.character(annot)))
})

if (requireNamespace('MuSiC', quietly = TRUE)) {
  test('music-hc leaves out the cell type of a single donor and returns the others', {
    groups <- celltype_group_map(as.character(annot), coarse)
    res <- suppressMessages(suppressWarnings(deconvolute_music_hc(bulk, counts, as.character(annot), batch, groups)))
    !('pDC_c4-IRF7+' %in% colnames(res)) && all(abs(rowSums(res) - 1) < 1e-6) && ncol(res) == length(unique(annot)) - 1
  })
}

if (requireNamespace('SECRET', quietly = TRUE)) {
  test('SECRET reference keeps cell types whose names contain each other apart', {
    ref <- suppressMessages(suppressWarnings(build_secret_reference(counts, as.character(annot), batch)))$reference
    # scRefer averages the cells of a cell type; with the names themselves 'Monocytes' would also
    # take the cells of 'Monocytes classical'
    genes <- rownames(ref)
    mono_only <- rowMeans(sweep(counts[genes, annot == 'Monocytes'], 2, colSums(counts[, annot == 'Monocytes']), '/'))
    all(c('Monocytes', 'Monocytes classical') %in% colnames(ref)) &&
      cor(ref[, 'Monocytes'], mono_only) > cor(ref[, 'Monocytes classical'], mono_only)
  })
}

test('split by donor and by cell: reproducible, invariants hold, factors become characters', {
  ok <- TRUE
  for (by in c('donor', 'cell')) {
    outs <- file.path(tmp_dir, paste0('split_', by, c('_a', '_b')))
    for (o in outs) {
      dir.create(o, showWarnings = FALSE)
      res <- system2('Rscript', c(file.path(pipeline_dir, 'bin', 'splitReferenceNF.R'), 'sc-ci', shQuote(dirname(sc_dir)),
                                  '0.5', '1', 'warn', '4', '1', shQuote(o), by), stdout = FALSE, stderr = FALSE)
      ok <- ok && res == 0
    }
    a <- read.csv(file.path(outs[1], 'split', 'assignment.csv'))
    ok <- ok && identical(a, read.csv(file.path(outs[2], 'split', 'assignment.csv')))
    donors_in_parts <- tapply(a$part, a$donor, function(p) length(unique(p)))
    ok <- ok && if (by == 'donor') all(donors_in_parts == 1) else any(donors_in_parts > 1)
    ok <- ok && is.character(readRDS(file.path(outs[1], 'deconvolution', 'sc-ci', 'celltype_annotations.rds')))
  }
  ok
})

if (length(problems) > 0) { cat(length(problems), 'unit test(s) failed\n'); quit(status = 1) }
cat('All unit tests passed\n')
