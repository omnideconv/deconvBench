#!/usr/bin/env Rscript
# Check that a trace file of the subsampling workflow can still be parsed the way
# visualisation/supplement/fig_s4.Rmd does it: tasks are selected by process name, and method,
# subset size and number of bulk samples are read from fixed word positions of the script column.
#
# Usage: check_trace_parsing.R <trace.csv>

library(data.table)

df <- fread(commandArgs(trailingOnly = TRUE)[1], fill = TRUE)

# same steps as in fig_s4.Rmd
df <- df[-which(df$task_id %in% c('\t', '')), ]
df <- df[rowSums(is.na(df)) != ncol(df), ]

script_lines <- seq(2, nrow(df), by = 2)
script_entries <- unlist(df[script_lines, 'task_id'])
df <- df[!script_lines]
df$script <- script_entries
df$script <- gsub('\'', '', df$script)
df$script <- gsub('\\t', ' ', df$script)
df$process <- gsub('subsampling:', '', df$process)

df_sig <- subset(df, process == 'CREATE_SIGNATURE_PREPROCESSED')
df_deconv <- subset(df, process == 'DECONVOLUTE')

df_deconv$method <- sapply(strsplit(df_deconv$script, ' '), `[[`, 7)
df_deconv$n_cells <- gsub('perc', '', sapply(strsplit(sapply(strsplit(df_deconv$script, ' '), `[[`, 3), '_'), `[[`, 2))
df_deconv$bulk_samples <- gsub('nsamples', '', sapply(strsplit(sapply(strsplit(df_deconv$script, ' '), `[[`, 5), '-'), `[[`, 3))

df_sig$method <- sapply(strsplit(df_sig$script, ' '), `[[`, 6)
df_sig$n_cells <- gsub('perc', '', sapply(strsplit(sapply(strsplit(df_sig$script, ' '), `[[`, 2), '_'), `[[`, 2))
df_sig$bulk_samples <- gsub('nsamples', '', sapply(strsplit(sapply(strsplit(df_sig$script, ' '), `[[`, 4), '-'), `[[`, 3))

parsed <- rbind(df_sig[, c('process', 'method', 'n_cells', 'bulk_samples')],
                df_deconv[, c('process', 'method', 'n_cells', 'bulk_samples')])
print(as.data.frame(parsed))

methods <- c('autogenes', 'bayesprism', 'bisque', 'cibersortx', 'dwls', 'music', 'rectangle', 'scaden', 'scdc')
ok <- nrow(df_sig) > 0 && nrow(df_deconv) > 0 &&
  all(parsed$method %in% methods) &&
  !anyNA(suppressWarnings(as.numeric(parsed$n_cells))) &&
  !anyNA(suppressWarnings(as.numeric(parsed$bulk_samples)))
if (!ok) {
  cat('Trace file cannot be parsed like fig_s4.Rmd expects.\n')
  quit(save = 'no', status = 1)
}
cat('Trace file is parsed as expected.\n')
