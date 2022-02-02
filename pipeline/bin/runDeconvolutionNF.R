#!/usr/bin/Rscript

"Usage: 
  runDeconvolutionNF.R <sc_path> <sc_datasetname> <rna_path> <rna_datasetname> <deconv_method> <signature> <remapping_sheet> [<coarse>]
Options:
<sc_path> path to sc dataset
<sc_datasetname> name of sc dataset
<rna_path> path to rnaseq dataset
<rna_datasetname> name of rnaseq dataset
<deconv_method>  deconv method
<signature> signature matrix
<remapping_sheet> excel mapping sheet with remapping
<coarse> logical, if TRUE celltypes are mapped to higher level" -> doc

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

signature <- readRDS(args$signature)

##remap celltype_annotations##
source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)


library(Biobase)#necessary for music
if(method=="cibersortx"){
  omnideconv::set_cibersortx_credentials("k.reinisch@campus.lmu.de", "a05114832330fda42ce0f5596875ee0d")
}

###autogenes ONLY works this way!! i don't get why###
# library("devtools")
# load_all("/nfs/proj/omnideconv_benchmarking/omnideconv")
# deconvolution <- deconvolute(bulk_gene_expression = get(paste(rnaseq_ds, "_pbmc_tpm", sep="")), 
#                              signature = signature, 
#                              single_cell_object = sc_matrix, 
#                              batch_ids = sc_batch, 
#                              cell_type_annotations = sc_celltype_annotations, 
#                              method = method)

runtime <- system.time({
deconvolution <- omnideconv::deconvolute(bulk_gene_expression = get(paste(rnaseq_ds, "_pbmc_tpm", sep="")), 
                        signature = signature, 
                        single_cell_object = sc_matrix, 
                        batch_ids = sc_batch, 
                        cell_type_annotations = sc_celltype_annotations, 
                        method = method)})
saveRDS(deconvolution, 
        paste("deconvolution_", method, "_", sc_ds, "_", rnaseq_ds, ".rds", sep=""), 
        compress = FALSE)
print(runtime)