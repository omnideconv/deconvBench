#!/usr/bin/env Rscript
 
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

sc_path <- args[2]
sc_ds <- args[4]
sc_matrix <- readRDS(file.path(sc_path, sc_ds, "matrix.rds"))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))

rnaseq_path <- args[6]
rnaseq_ds <- args[8]
load(file.path(rnaseq_path, rnaseq_ds, paste(rnaseq_ds, "_pbmc_tpm.RData", sep="")))

method <- args[10]

signature <- readRDS(args[12])

library(Biobase)#necessary for music
if(method=="cibersortx"){
  omnideconv::set_cibersortx_credentials("k.reinisch@campus.lmu.de", "a05114832330fda42ce0f5596875ee0d")
}


deconvolution <- omnideconv::deconvolute(bulk_gene_expression = get(paste(rnaseq_ds, "_pbmc_tpm", sep="")), 
                        signature = signature, 
                        single_cell_object = sc_matrix, 
                        batch_ids = sc_batch, 
                        cell_type_annotations = sc_celltype_annotations, 
                        method = method)

saveRDS(deconvolution, 
        paste("deconvolution_", method, "_", sc_ds, "_", rnaseq_ds, ".rds", sep=""), 
        compress = FALSE)
