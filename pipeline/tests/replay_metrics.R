#!/usr/bin/env Rscript
# Run two versions of computeMetricsNF.R on stored deconvolution results and compare what they
# compute. This covers the dataset-specific cell-type sums (hoek, morandini, altman, dietrich),
# which the fixture does not trigger. Has to run inside the benchmark container.
#
# Usage: replay_metrics.R <results_dir> <old_pipeline_dir> <new_pipeline_dir> <scratch_dir> <bulk_base_dir>...
#
#   <results_dir>        stored results of the main workflow, one folder per result
#   <old_pipeline_dir>   pipeline directory whose script builds the result folder name itself
#   <new_pipeline_dir>   pipeline directory whose script is given the result folder
#   <bulk_base_dir>...   directories that contain <bulk_name>/<bulk_name>_facs.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: replay_metrics.R <results_dir> <old_pipeline_dir> <new_pipeline_dir> <scratch_dir> <bulk_base_dir>...")
results_dir <- args[1]
old_dir <- args[2]
new_dir <- args[3]
scratch <- args[4]
bulk_bases <- args[-(1:4)]

needed <- c("deconvolution.rds", "runtime_signature.rds", "runtime_deconvolution.rds")
status <- c()

for (folder in list.dirs(results_dir, recursive = FALSE, full.names = FALSE)) {
  src <- file.path(results_dir, folder)
  if (!all(file.exists(file.path(src, c(needed, "results_metric.rds"))))) next

  # method and dataset names of this result are stored in its runtimes table
  runtimes <- readRDS(file.path(src, "results_metric.rds"))$runtimes
  method <- as.character(runtimes$method[1])
  sc_name <- as.character(runtimes$sc_ds[1])
  bulk_name <- as.character(runtimes$bulk_ds[1])
  bulk_base <- bulk_bases[file.exists(file.path(bulk_bases, bulk_name, paste0(bulk_name, "_facs.rds")))][1]
  if (is.na(bulk_base)) {
    status[folder] <- "skipped: ground truth not found"
    next
  }

  old_res <- file.path(scratch, "old", folder)
  new_res <- file.path(scratch, "new", folder)
  for (d in c(old_res, new_res)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(src, needed), d, overwrite = TRUE)
  }

  run <- function(pipeline_dir, result_arg) {
    system2(file.path(pipeline_dir, "bin", "computeMetricsNF.R"),
            shQuote(c(sc_name, "sc", bulk_name, bulk_base, method, "0", "0", result_arg, pipeline_dir)),
            stdout = FALSE, stderr = FALSE)
  }
  exit_old <- run(old_dir, file.path(scratch, "old"))
  exit_new <- run(new_dir, new_res)

  if (exit_old != 0 || exit_new != 0) {
    status[folder] <- if (exit_old == exit_new) "both versions fail" else paste0("DIFFERENT: exit status old ", exit_old, ", new ", exit_new)
    next
  }
  old <- readRDS(file.path(old_res, "results_metric.rds"))
  new <- readRDS(file.path(new_res, "results_metric.rds"))
  # all.equal with tolerance 0 demands the same numbers but, unlike identical(), does not compare
  # the address of the pandas index that results of Python methods carry as an attribute
  comparison <- all.equal(old, new, tolerance = 0)
  status[folder] <- if (identical(old, new) || isTRUE(comparison)) "identical" else paste("DIFFERENT:", paste(utils::head(comparison, 2), collapse = " | "))
}

bulk_of <- sub(".*_(counts|cpm)_(.*)_(counts|tpm)_ct.*", "\\2", names(status))
print(table(bulk = bulk_of, status = sub(":.*", "", status)))
bad <- status[grepl("^DIFFERENT", status)]
if (length(bad) > 0) {
  for (n in names(bad)) cat(n, "\n   ", bad[n], "\n")
  quit(save = "no", status = 1)
}
cat("\nBoth versions compute the same metrics for", sum(status == "identical"), "results.\n")
