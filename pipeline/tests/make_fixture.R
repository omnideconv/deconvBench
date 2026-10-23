#!/usr/bin/env Rscript
# Build the regression-test fixture from the small datasets tracked in data/.
# Has to run inside the benchmark container (needs SimBu).
#
# Usage: make_fixture.R <repo_data_dir> <fixture_dir>
#
# Result:
#   <fixture_dir>/single_cell/sc1    all files needed by every workflow (simbu_ds, fine/normal/coarse annotations)
#   <fixture_dir>/single_cell/sc1nc  same as sc1 without coarse annotation (resolution workflow, no-coarse branch)
#   <fixture_dir>/single_cell/sc2    second reference with fewer cells and cell types (impact technology workflow)
#   <fixture_dir>/bulks/bulk1|bulk2  copies of the tracked bulks
#   <fixture_dir>/bulks/hoek         copy of bulk1 under the name that triggers the hoek-specific code

library(SimBu)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: make_fixture.R <repo_data_dir> <fixture_dir>")
data_dir <- args[1]
fixture_dir <- args[2]

sc_in <- file.path(data_dir, "single_cell", "sc1")
sc_out <- file.path(fixture_dir, "single_cell")
bulk_out <- file.path(fixture_dir, "bulks")
dir.create(sc_out, recursive = TRUE, showWarnings = FALSE)
dir.create(bulk_out, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(sc_in, "matrix_counts.rds"))
norm <- readRDS(file.path(sc_in, "matrix_norm_counts.rds"))
annot <- readRDS(file.path(sc_in, "celltype_annotations.rds"))
coarse <- readRDS(file.path(sc_in, "celltype_annotations_coarse.rds"))
batch <- readRDS(file.path(sc_in, "batch.rds"))

write_sc <- function(name, idx, with_coarse = TRUE, with_resolution = TRUE) {
  out <- file.path(sc_out, name)
  dir.create(out, showWarnings = FALSE)
  cur_counts <- counts[, idx]
  cur_norm <- norm[, idx]
  cur_annot <- annot[idx]
  saveRDS(cur_counts, file.path(out, "matrix_counts.rds"))
  saveRDS(cur_norm, file.path(out, "matrix_norm_counts.rds"))
  saveRDS(cur_annot, file.path(out, "celltype_annotations.rds"))
  saveRDS(batch[idx], file.path(out, "batch.rds"))

  simbu_ds <- SimBu::dataset(
    annotation = data.frame(ID = colnames(cur_counts), cell_type = cur_annot),
    count_matrix = cur_counts,
    tpm_matrix = cur_norm,
    name = name
  )
  saveRDS(simbu_ds, file.path(out, "simbu_ds.rds"))

  if (with_resolution) {
    # fine annotation: split three cell types into two subtypes each, alternating by position
    fine <- cur_annot
    split_types <- c("B cells" = "B cells", "Monocytes" = "Monocytes", "T cells CD8" = "T cells CD8")
    for (ct in names(split_types)) {
      pos <- which(cur_annot == ct)
      fine[pos] <- paste0(ct, c(" s1", " s2"))[(seq_along(pos) %% 2) + 1]
    }
    saveRDS(cur_annot, file.path(out, "celltype_annotations_normal.rds"))
    saveRDS(fine, file.path(out, "celltype_annotations_fine.rds"))
    mapping <- unique(data.frame(Fine = fine, Normal = cur_annot, Coarse = coarse[idx]))
    mapping <- mapping[order(mapping$Fine), ]
    write.table(mapping, file.path(out, "cell_type_mappings.csv"), sep = ",", row.names = FALSE, quote = FALSE)
    if (with_coarse) {
      saveRDS(coarse[idx], file.path(out, "celltype_annotations_coarse.rds"))
    }
  }
}

all_cells <- seq_along(annot)
write_sc("sc1", all_cells)
write_sc("sc1nc", all_cells, with_coarse = FALSE)
# second reference: two thirds of the cells, without two of the cell types
write_sc("sc2", which(all_cells %% 3 != 0 & !(annot %in% c("Platelet", "ILC"))), with_resolution = FALSE)

for (b in c("bulk1", "bulk2")) {
  dir.create(file.path(bulk_out, b), showWarnings = FALSE)
  file.copy(list.files(file.path(data_dir, "bulks", b), full.names = TRUE), file.path(bulk_out, b), overwrite = TRUE)
}
dir.create(file.path(bulk_out, "hoek"), showWarnings = FALSE)
for (suffix in c("counts", "tpm", "facs")) {
  file.copy(file.path(data_dir, "bulks", "bulk1", paste0("bulk1_", suffix, ".rds")),
            file.path(bulk_out, "hoek", paste0("hoek_", suffix, ".rds")), overwrite = TRUE)
}

cat("Fixture written to", fixture_dir, "\n")
