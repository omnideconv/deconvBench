# Check cell_type_mappings.csv of single-cell datasets against their annotation files: every fine
# cell type has to be in the mapping, with the normal and coarse cell type of its cells. The
# resolution analysis sums the fine fractions with this table.
#
# Usage: Rscript check_mappings.R <dataset_dir> [<dataset_dir> ...]   (exit status 1 on a problem)

check_mapping <- function(dir) {
  f <- file.path(dir, 'cell_type_mappings.csv')
  if (!file.exists(f) || !file.exists(file.path(dir, 'celltype_annotations_fine.rds'))) return(character())
  m <- read.csv(f, check.names = FALSE)
  levels <- c(Fine = 'fine', Normal = 'normal', Coarse = 'coarse')
  levels <- levels[file.exists(file.path(dir, paste0('celltype_annotations_', levels, '.rds')))]
  observed <- unique(as.data.frame(lapply(setNames(levels, names(levels)), function(l)
    as.character(readRDS(file.path(dir, paste0('celltype_annotations_', l, '.rds')))))))
  problems <- character()
  missing <- setdiff(observed$Fine, m$Fine)
  if (length(missing) > 0) problems <- c(problems, paste0(basename(dir), ': fine cell types not in cell_type_mappings.csv: ', paste(missing, collapse = ', ')))
  extra <- setdiff(m$Fine, observed$Fine)
  if (length(extra) > 0) problems <- c(problems, paste0(basename(dir), ': cell_type_mappings.csv lists fine cell types that no cell has: ', paste(extra, collapse = ', ')))
  joined <- merge(observed, m, by = 'Fine', suffixes = c('', '.mapping'))
  for (l in setdiff(names(levels), 'Fine')) {
    wrong <- joined$Fine[joined[[l]] != joined[[paste0(l, '.mapping')]]]
    if (length(wrong) > 0) problems <- c(problems, paste0(basename(dir), ': ', l, ' level differs from the annotation for: ', paste(wrong, collapse = ', ')))
  }
  problems
}

if (sys.nframe() == 0) {
  problems <- unlist(lapply(commandArgs(trailingOnly = TRUE), check_mapping))
  if (length(problems) > 0) { cat(problems, sep = '\n'); quit(status = 1) }
  cat('mappings ok\n')
}
