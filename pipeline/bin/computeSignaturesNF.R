#!/usr/bin/env Rscript
'
computeSignatureNF.R
' 
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
#args <- list("","/nfs/data/omnideconv_benchmarking/data/singleCell", "","lambrechts","", "/nfs/data/omnideconv_benchmarking/data/PBMC","","hoek","","bseqsc")

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


if(method=="cibersortx"){
  omnideconv::set_cibersortx_credentials("k.reinisch@campus.lmu.de", "a05114832330fda42ce0f5596875ee0d")
}

signature <- omnideconv::build_model(single_cell_object = as.data.frame(sc_matrix), 
                                     cell_type_annotations = sc_celltype_annotations, 
				                            bulk_gene_expression = get(paste(rnaseq_ds, "_pbmc_tpm", sep="")), 
				                            markers = sc_marker,
				                            batch_ids = sc_batch, 
				                            verbose=TRUE,
				                            method = method)
print(signature)
saveRDS(signature, 
        paste("signature_", method, "_", sc_ds, "_", rnaseq_ds, ".rds", sep=""), 
        compress = FALSE)
