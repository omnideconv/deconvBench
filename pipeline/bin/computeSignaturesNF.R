#!/usr/bin/Rscript
print("start")
library(conflicted)
conflicted::conflict_scout()
#library(docopt)

"Usage: 
  computeSignatureNF.R <sc_path> <sc_datasetname> <rna_path> <rna_datasetname> <deconv_method> <remapping_sheet>
Options:
<sc_path> path to sc dataset
<sc_datasetname> name of sc dataset
<rna_path> path to rnaseq dataset
<rna_datasetname> name of rnaseq dataset
<deconv_method>  deconv method
<remapping_sheet> excel mapping sheet with remapping" -> doc

args <- docopt::docopt(doc)

sc_path <- args$sc_path
sc_ds <- args$sc_datasetname
sc_matrix <- readRDS(file.path(sc_path, sc_ds, "matrix.rds"))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))

rnaseq_path <- args$rna_path
rnaseq_ds <- args$rna_datasetname
load(file.path(rnaseq_path, rnaseq_ds, paste(rnaseq_ds, "_pbmc_tpm.RData", sep="")))

method <- args$deconv_method
remapping_sheet <- args$remapping_sheet

source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)

if(method=="cibersortx"){
  omnideconv::set_cibersortx_credentials("k.reinisch@campus.lmu.de", 
                                         "a05114832330fda42ce0f5596875ee0d")
}

runtime <- system.time({
signature <- omnideconv::build_model(single_cell_object = as.data.frame(sc_matrix), 
                                    cell_type_annotations = sc_celltype_annotations, 
				                            bulk_gene_expression = get(paste(rnaseq_ds, "_pbmc_tpm", sep="")), 
				                            markers = sc_marker,
				                            batch_ids = sc_batch, 
				                            verbose=TRUE,
				                            method = method)})

if(method %in% c("autogenes", "scaden")){ 
  #copy signature into working dir since Rtemp dir will be closed
  file.copy(signature, getwd(), recursive = TRUE)
 #rename signature
 signature <- file.path(getwd(), basename(signature))
}
saveRDS(signature, 
       paste("signature_", method, "_", sc_ds, "_", rnaseq_ds, ".rds", sep=""), 
        compress = FALSE)
print(runtime)