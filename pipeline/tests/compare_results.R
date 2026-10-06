#!/usr/bin/env Rscript
# Compare two result trees of the pipeline (e.g. reference run vs. run of refactored code).
#
# Usage: compare_results.R <reference_dir> <new_dir> [--tolerance=<x>] [--report=<file.csv>]
#
# 1) both trees have to contain the same files (relative paths)
# 2) every .rds/.RData/.csv has to have equal content
#
# Not compared: tmp/ directories (method scratch), runtime_*.rds (timings) and the content of
# scaden model/ directories (only their presence is checked).
# Exit status 0 if the trees are equivalent, 1 otherwise.

args <- commandArgs(trailingOnly = TRUE)
opts <- args[grepl("^--", args)]
dirs <- args[!grepl("^--", args)]
if (length(dirs) != 2) stop("Usage: compare_results.R <reference_dir> <new_dir> [--tolerance=<x>] [--report=<file.csv>]")

get_opt <- function(name, default) {
  hit <- opts[startsWith(opts, paste0("--", name, "="))]
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
tolerance <- as.numeric(get_opt("tolerance", "0"))
report_file <- get_opt("report", NA)

ref_dir <- normalizePath(dirs[1])
new_dir <- normalizePath(dirs[2])

list_tree <- function(d) {
  f <- list.files(d, recursive = TRUE, all.files = TRUE)
  f <- f[!grepl("(^|/)tmp/", f)]
  f <- f[!grepl("(^|/)tmp_(input|output)/", f)]
  sort(f)
}

content_checked <- function(f) {
  !grepl("(^|/)runtime_[^/]*\\.rds$", f) & !grepl("(^|/)model/", f)
}

compare_objects <- function(a, b) {
  if (identical(a, b)) return("identical")
  res <- tryCatch(
    if (tolerance > 0) all.equal(a, b, tolerance = tolerance) else all.equal(a, b, tolerance = 0),
    error = function(e) paste("all.equal failed:", conditionMessage(e))
  )
  if (isTRUE(res)) {
    if (tolerance > 0) return("equal within tolerance")
    return("identical")
  }
  paste("DIFFERENT:", paste(utils::head(res, 3), collapse = " | "))
}

compare_file <- function(f) {
  a <- file.path(ref_dir, f)
  b <- file.path(new_dir, f)
  if (!content_checked(f)) return("not compared")
  if (grepl("\\.rds$", f, ignore.case = TRUE)) {
    return(compare_objects(readRDS(a), readRDS(b)))
  }
  if (grepl("\\.RData$", f, ignore.case = TRUE)) {
    ea <- new.env(); eb <- new.env()
    load(a, envir = ea); load(b, envir = eb)
    return(compare_objects(as.list(ea)[sort(ls(ea))], as.list(eb)[sort(ls(eb))]))
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

if (!is.na(report_file)) write.csv(report, report_file, row.names = FALSE)

cat("reference:", ref_dir, "\nnew:      ", new_dir, "\ntolerance:", tolerance, "\n\n")
summary_tbl <- table(report$workflow, ifelse(report$ok, sub(":.*", "", report$status), sub(":.*", "", report$status)))
print(summary_tbl)

bad <- report[!report$ok, ]
if (nrow(bad) > 0) {
  cat("\n", nrow(bad), " file(s) differ:\n", sep = "")
  for (i in seq_len(nrow(bad))) cat("  ", bad$file[i], "\n      ", bad$status[i], "\n", sep = "")
  quit(save = "no", status = 1)
}
cat("\nResult trees are equivalent (", nrow(report), " files).\n", sep = "")
