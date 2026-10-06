#!/usr/bin/env Rscript
# Compare two result trees of the pipeline (e.g. reference run vs. run of refactored code).
#
# Usage: compare_results.R <reference_dir> <new_dir> [--tolerance=<x>] [--report=<file.csv>] [--known-unstable=<file>] [--only=<regex>]
#
# 1) both trees have to contain the same files (relative paths)
# 2) every .rds/.RData/.csv has to have equal content, up to a numeric tolerance (default 1e-8;
#    two runs of the same code on different cluster nodes differ by about 1e-13)
#
# Not compared: tmp/ directories (method scratch), runtime_*.rds and the runtimes entry of
# results_metric.rds (timings), and the content of scaden model/ directories (presence only).
# Scaden is not reproducible (no seed can be passed through omnideconv), so for scaden results
# only the structure is compared: dimensions, row names and the set of column names.
# Files matching a pattern in the known-unstable list (default: known_unstable.txt next to this
# script) are reported but do not fail the comparison.
# With --only, the comparison is restricted to files whose relative path matches the regular
# expression, e.g. --only='^(main|simulation)/' to check some workflows only.
# Exit status 0 if the trees are equivalent, 1 otherwise.

args <- commandArgs(trailingOnly = TRUE)
opts <- args[grepl("^--", args)]
dirs <- args[!grepl("^--", args)]
if (length(dirs) != 2) stop("Usage: compare_results.R <reference_dir> <new_dir> [--tolerance=<x>] [--report=<file.csv>]")

get_opt <- function(name, default) {
  hit <- opts[startsWith(opts, paste0("--", name, "="))]
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
tolerance <- as.numeric(get_opt("tolerance", "1e-8"))
report_file <- get_opt("report", NA)
only <- get_opt("only", NA)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
unstable_file <- get_opt("known-unstable", file.path(dirname(script_path), "known_unstable.txt"))
unstable_patterns <- character(0)
if (file.exists(unstable_file)) {
  unstable_patterns <- readLines(unstable_file)
  unstable_patterns <- unstable_patterns[nzchar(trimws(unstable_patterns)) & !startsWith(unstable_patterns, "#")]
}

ref_dir <- normalizePath(dirs[1])
new_dir <- normalizePath(dirs[2])

list_tree <- function(d) {
  f <- list.files(d, recursive = TRUE, all.files = TRUE)
  f <- f[!grepl("(^|/)tmp/", f)]
  f <- f[!grepl("(^|/)tmp_(input|output)/", f)]
  if (!is.na(only)) f <- f[grepl(only, f, perl = TRUE)]
  sort(f)
}

content_checked <- function(f) {
  !grepl("(^|/)runtime_[^/]*\\.rds$", f) & !grepl("(^|/)model/", f)
}

is_scaden <- function(f) {
  grepl("(^|/)scaden(_|/)", f)
}

# reduce an object to what has to be stable for a stochastic method
structure_of <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(list(dim = dim(x), rows = rownames(x), cols = sort(colnames(x))))
  }
  if (is.list(x)) {
    return(lapply(x[sort(names(x))], structure_of))
  }
  list(class = class(x), length = length(x))
}

compare_objects <- function(a, b, structure_only = FALSE) {
  if (identical(a, b)) return("identical")
  if (structure_only) {
    if (identical(structure_of(a), structure_of(b))) return("same structure")
    return("DIFFERENT: structure")
  }
  res <- tryCatch(
    all.equal(a, b, tolerance = tolerance),
    error = function(e) paste("all.equal failed:", conditionMessage(e))
  )
  if (isTRUE(res)) return("equal within tolerance")
  paste("DIFFERENT:", paste(utils::head(res, 3), collapse = " | "))
}

compare_file <- function(f) {
  a <- file.path(ref_dir, f)
  b <- file.path(new_dir, f)
  if (!content_checked(f)) return("not compared")
  if (grepl("\\.rds$", f, ignore.case = TRUE)) {
    obj_a <- readRDS(a)
    obj_b <- readRDS(b)
    if (grepl("(^|/)results_metric\\.rds$", f)) {
      obj_a$runtimes <- NULL
      obj_b$runtimes <- NULL
    }
    return(compare_objects(obj_a, obj_b, structure_only = is_scaden(f)))
  }
  if (grepl("\\.RData$", f, ignore.case = TRUE)) {
    ea <- new.env(); eb <- new.env()
    load(a, envir = ea); load(b, envir = eb)
    return(compare_objects(as.list(ea)[sort(ls(ea))], as.list(eb)[sort(ls(eb))], structure_only = is_scaden(f)))
  }
  if (grepl("\\.csv$", f, ignore.case = TRUE)) {
    return(compare_objects(read.csv(a), read.csv(b), structure_only = is_scaden(f)))
  }
  if (unname(tools::md5sum(a)) == unname(tools::md5sum(b))) "identical" else "DIFFERENT: file content"
}

ref_files <- list_tree(ref_dir)
new_files <- list_tree(new_dir)

only_ref <- setdiff(ref_files, new_files)
only_new <- setdiff(new_files, ref_files)
common <- intersect(ref_files, new_files)

status <- vapply(common, compare_file, character(1))
report <- rbind(
  data.frame(file = common, status = unname(status)),
  data.frame(file = only_ref, status = rep("MISSING in new", length(only_ref))),
  data.frame(file = only_new, status = rep("ADDITIONAL in new", length(only_new)))
)
report$workflow <- sub("/.*", "", report$file)
report$ok <- !grepl("^(DIFFERENT|MISSING|ADDITIONAL)", report$status)
unstable <- Reduce(`|`, lapply(unstable_patterns, grepl, x = report$file), rep(FALSE, nrow(report)))
report$status[unstable & !report$ok] <- paste("known unstable,", report$status[unstable & !report$ok])
report$ok <- report$ok | unstable

if (!is.na(report_file)) write.csv(report, report_file, row.names = FALSE)

cat("reference:", ref_dir, "\nnew:      ", new_dir, "\ntolerance:", tolerance, "\n\n")
summary_tbl <- table(report$workflow, sub("[:,].*", "", report$status))
print(summary_tbl)

bad <- report[!report$ok, ]
if (nrow(bad) > 0) {
  cat("\n", nrow(bad), " file(s) differ:\n", sep = "")
  for (i in seq_len(nrow(bad))) cat("  ", bad$file[i], "\n      ", bad$status[i], "\n", sep = "")
  quit(save = "no", status = 1)
}
cat("\nResult trees are equivalent (", nrow(report), " files).\n", sep = "")
