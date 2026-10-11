#!/usr/bin/Rscript

print("Starting reference split script ...")

library(docopt)
library(Matrix)

"Usage:
  splitReferenceNF.R <sc_name> <sc_path> <fraction_simulation> <min_cells> <min_cells_action> <seed> <reference_cells_fraction> <output_dir>
Options:
<sc_name> name of sc dataset
<sc_path> path to sc dataset
<fraction_simulation> share of the donors whose cells are used for the simulation of pseudo-bulks
<min_cells> minimum number of cells that every cell type needs in each part
<min_cells_action> what to do if a cell type has fewer cells in a part: error or warn
<seed> seed of the donor assignment
<reference_cells_fraction> share of the cells of the method part that is kept (1 keeps all)
<output_dir> output directory" -> doc

args <- docopt::docopt(doc)

sc_dataset <- args$sc_name
sc_dir <- file.path(args$sc_path, sc_dataset)
fraction_simulation <- as.numeric(args$fraction_simulation)
min_cells <- as.numeric(args$min_cells)
min_cells_action <- args$min_cells_action
seed <- as.integer(args$seed)
reference_cells_fraction <- as.numeric(args$reference_cells_fraction)
output_dir <- args$output_dir

if(is.na(fraction_simulation) || fraction_simulation <= 0 || fraction_simulation >= 1){
  stop('split_fraction_simulation has to be between 0 and 1')
}
if(is.na(reference_cells_fraction) || reference_cells_fraction <= 0 || reference_cells_fraction > 1){
  stop('split_reference_cells_fraction has to be above 0 and at most 1')
}
if(!min_cells_action %in% c('error', 'warn')){
  stop('split_min_cells_action has to be error or warn')
}

annotations <- as.character(readRDS(file.path(sc_dir, 'celltype_annotations.rds')))
n_cells <- length(annotations)

# The donor of every cell: donor.rds if the dataset has one, otherwise batch.rds. A dataset whose
# batches are samples of the same donors needs a donor.rds, or cells of one donor end up in both parts.
donor_file <- if(file.exists(file.path(sc_dir, 'donor.rds'))) 'donor.rds' else 'batch.rds'
donors <- as.character(readRDS(file.path(sc_dir, donor_file)))
if(length(donors) != n_cells){
  stop(paste0(donor_file, ' has ', length(donors), ' entries, celltype_annotations.rds has ', n_cells))
}
donor_ids <- sort(unique(donors))
if(length(donor_ids) < 2){
  stop(paste0('The dataset has only one donor in ', donor_file, ' and cannot be split'))
}

###########################
#### assign the donors ####
###########################

set.seed(seed)
shuffled <- sample(donor_ids)
n_simulation <- min(max(round(fraction_simulation * length(donor_ids)), 1), length(donor_ids) - 1)
simulation_donors <- shuffled[seq_len(n_simulation)]
part <- ifelse(donors %in% simulation_donors, 'simulation', 'deconvolution')

# optionally keep only a share of the cells of the method part, drawn within every cell type
keep <- rep(TRUE, n_cells)
if(reference_cells_fraction < 1){
  for(cell_type in sort(unique(annotations))){
    candidates <- which(part == 'deconvolution' & annotations == cell_type)
    if(length(candidates) > 1){
      n_keep <- max(round(reference_cells_fraction * length(candidates)), 1)
      keep[setdiff(candidates, sample(candidates, n_keep))] <- FALSE
    }
  }
}
part[!keep] <- 'not used'

assignment <- data.frame(cell = seq_len(n_cells), donor = donors, cell_type = annotations, part = part, stringsAsFactors = FALSE)
summary_table <- as.data.frame.matrix(table(factor(annotations), factor(part, levels = c('simulation', 'deconvolution', 'not used'))))
summary_table <- cbind(cell_type = rownames(summary_table), summary_table)
donor_table <- data.frame(donor = donor_ids, part = ifelse(donor_ids %in% simulation_donors, 'simulation', 'deconvolution'),
                          cells = as.vector(table(factor(donors, levels = donor_ids))))

cat(paste0('Donors (', donor_file, '): ', n_simulation, ' for the simulation, ', length(donor_ids) - n_simulation, ' for the methods\n'))
print(summary_table, row.names = FALSE)

too_few <- summary_table$cell_type[summary_table$simulation < min_cells | summary_table$deconvolution < min_cells]
if(length(too_few) > 0){
  message_text <- paste0('Fewer than ', min_cells, ' cells in one part of the split of ', sc_dataset, ' for: ', paste(too_few, collapse = ', '),
                         '. Change split_seed, split_fraction_simulation or split_min_cells, or set split_min_cells_action to warn.')
  if(min_cells_action == 'error'){
    stop(message_text)
  }
  warning(message_text, immediate. = TRUE)
}

##########################
#### write both parts ####
##########################

split_dir <- file.path(output_dir, 'split')
dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(assignment, file.path(split_dir, 'assignment.csv'), row.names = FALSE)
write.csv(summary_table, file.path(split_dir, 'cells_per_cell_type.csv'), row.names = FALSE)
write.csv(donor_table, file.path(split_dir, 'donors.csv'), row.names = FALSE)
write.csv(data.frame(parameter = c('dataset', 'donors_from', 'split_fraction_simulation', 'split_min_cells', 'split_min_cells_action', 'split_seed', 'split_reference_cells_fraction'),
                     value = c(sc_dataset, donor_file, fraction_simulation, min_cells, min_cells_action, seed, reference_cells_fraction)),
          file.path(split_dir, 'parameters.csv'), row.names = FALSE)

part_dirs <- list(simulation = file.path(output_dir, 'simulation', sc_dataset), deconvolution = file.path(output_dir, 'deconvolution', sc_dataset))
for(part_dir in part_dirs){
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
}

# one value per cell: annotations of every level, batch and donor
per_cell_files <- c(list.files(sc_dir, pattern = '^celltype_annotations.*\\.rds$'), 'batch.rds', if(donor_file == 'donor.rds') 'donor.rds')
for(file_name in per_cell_files){
  values <- readRDS(file.path(sc_dir, file_name))
  if(length(values) != n_cells){
    cat(paste0('Not split, because it does not have one entry per cell: ', file_name, '\n'))
    next
  }
  for(part_name in names(part_dirs)){
    saveRDS(values[part == part_name], file.path(part_dirs[[part_name]], file_name))
  }
}

# expression matrices, one column per cell
for(file_name in c('matrix_counts.rds', 'matrix_norm_counts.rds')){
  expression <- readRDS(file.path(sc_dir, file_name))
  if(ncol(expression) != n_cells){
    stop(paste0(file_name, ' has ', ncol(expression), ' columns, celltype_annotations.rds has ', n_cells, ' entries'))
  }
  for(part_name in names(part_dirs)){
    saveRDS(expression[, part == part_name], file.path(part_dirs[[part_name]], file_name))
  }
  rm(expression)
  gc()
}

# Prepared SimBu dataset of the simulation part. Its cells are in the order of the expression matrices.
simbu_file <- file.path(sc_dir, 'simbu_ds.rds')
if(file.exists(simbu_file)){
  simbu_ds <- readRDS(simbu_file)
  if(ncol(simbu_ds) != n_cells){
    stop(paste0('simbu_ds.rds has ', ncol(simbu_ds), ' cells, celltype_annotations.rds has ', n_cells, ' entries'))
  }
  if(!all(as.character(SummarizedExperiment::colData(simbu_ds)$cell_type) == annotations)){
    stop('The cell types in simbu_ds.rds are not in the order of celltype_annotations.rds')
  }
  saveRDS(simbu_ds[, part == 'simulation'], file.path(part_dirs$simulation, 'simbu_ds.rds'))
}

# files that do not depend on the cells
for(file_name in intersect('cell_type_mappings.csv', list.files(sc_dir))){
  for(part_dir in part_dirs){
    file.copy(file.path(sc_dir, file_name), file.path(part_dir, file_name))
  }
}

cat(paste0("Split of '", sc_dataset, "' stored in '", output_dir, "'\n"))
