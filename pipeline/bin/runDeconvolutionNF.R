#!/usr/bin/Rscript

"Usage: 
  runDeconvolutionNF.R <sc_matrix> <sc_meta> <sc_path> <rna_path> <rna_datasetname> <deconv_method> <signature> <remapping_sheet> [<coarse>]
Options:
<sc_matrix> path to sc matrix
<sc_meta> meta information of sc dataset
<sc_path> path to sc dataset
<rna_path> path to rnaseq dataset
<rna_datasetname> name of rnaseq dataset
<deconv_method>  deconv method
<signature> signature matrix
<remapping_sheet> excel mapping sheet with remapping
<coarse> logical, if TRUE celltypes are mapped to higher level" -> doc

args <- docopt::docopt(doc)
sc_path <- args$sc_path
sc_meta <- strsplit(gsub("\\]", "", gsub("\\[", "", args$sc_meta)), split = ", ")[[1]]
sc_ds <- strsplit(sc_meta[1], split = ":")[[1]][2]
sc_type <- strsplit(sc_meta[2], ":")[[1]][2]


sc_matrix <- readRDS(file.path(args$sc_matrix))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))

rnaseq_path <- args$rna_path
rnaseq_ds <- args$rna_datasetname
rnaseq_data <- readRDS(file.path(rnaseq_path, rnaseq_ds, paste(rnaseq_ds, "_", sc_type, ".rds", sep="")))

method <- args$deconv_method
remapping_sheet <- args$remapping_sheet

signature <- readRDS(args$signature)

##remap celltype_annotations##
source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)

escapeCelltypesAutogenes <- function(celltype){
  return(gsub(" ", "xxxx", celltype))
}
reEscapeCelltypesAutogenes <- function(celltype){
  celltype <- gsub("21b2c6e87f8711ec9bf265fb9bf6ab9c", "+", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9a", "-", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9f", "(", celltype)
  celltype <- gsub("21b2c7567f8711ec9bf265fb9bf6ab9g", ")", celltype)
  return(gsub("xxxx", " ", celltype))
}


library(Biobase)#necessary for music

###autogenes ONLY works this way!! i don't get why###
if(method=="autogenes"){
  library("devtools")
  load_all("/nfs/proj/omnideconv_benchmarking/omnideconv")
}


if(method=="cibersortx"){
  omnideconv::set_cibersortx_credentials("k.reinisch@campus.lmu.de", "a05114832330fda42ce0f5596875ee0d")
} else if(method =="scaden" & sc_type=="counts"){
  saveRDS(NULL, paste("deconvolution_", method, "_", sc_ds, "_", sc_type, "_", rnaseq_ds, ".rds", sep=""), 
         compress = FALSE)
  quit()
} else if(method=="autogenes"){
  sc_celltype_annotations <- escapeCelltypesAutogenes(sc_celltype_annotations)
}



runtime <- system.time({
if(method=="cpm"){
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = rnaseq_data, 
                                           signature = signature, 
                                           single_cell_object = sc_matrix, 
                                           batch_ids = sc_batch, 
                                           cell_type_annotations = sc_celltype_annotations, 
                                           method = method)
                                           #no_cores = 120)
} else {
  deconvolution <- omnideconv::deconvolute(bulk_gene_expression = rnaseq_data, 
                        signature = signature, 
                        single_cell_object = sc_matrix, 
                        batch_ids = sc_batch, 
                        cell_type_annotations = sc_celltype_annotations, 
                        method = method)
  if(method=="autogenes"){
    colnames(deconvolution) <- reEscapeCelltypesAutogenes(colnames(deconvolution))
  }
}
})
saveRDS(deconvolution, 
        paste("deconvolution_", method, "_", sc_ds, "_", sc_type, "_", rnaseq_ds, ".rds", sep=""), 
        compress = FALSE)
print(runtime)
