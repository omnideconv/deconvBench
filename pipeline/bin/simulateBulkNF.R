#!/usr/bin/Rscript
print("start simulator")

"Usage: 
  simulateBulkNF.R <sc_path> <sc_datasetname> <rna_path> <rna_datasetname> <deconv_method> <remapping_sheet>
Options:
<sc_path> path to sc dataset
<sc_datasetname> name of sc dataset
<rna_path> path to rnaseq dataset
<rna_datasetname> name of rnaseq dataset
<deconv_method>  deconv method
<remapping_sheet> excel mapping sheet with remapping" -> doc

print(doc)
args <- docopt::docopt(doc)

sc_path <- args$sc_path
sc_ds <- args$sc_datasetname
sc_matrix <- readRDS(file.path(sc_path, sc_ds, "matrix.rds"))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))

method <- args$deconv_method
remapping_sheet <- args$remapping_sheet

source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)

library(SimBu)
###uniform###
simulated_bulk <- simulate_bulk(dataset(annotation = data.frame(ID = colnames(sc_matrix), cell_type = sc_celltype_annotations), 
                                                count_matrix = sc_matrix,
                                                name = sc_ds))
signature <- omnideconv::build_model(single_cell_object = as.data.frame(sc_matrix), 
                                     cell_type_annotations = sc_celltype_annotations, 
                                     bulk_gene_expression = simulated_bulk$pseudo_bulk, 
                                     markers = sc_marker,
                                     batch_ids = sc_batch, 
                                     verbose=TRUE,
                                     method = method)

deconvolution <- omnideconv::deconvolute(bulk_gene_expression = simulated_bulk$pseudo_bulk, 
                                         signature = signature, 
                                         single_cell_object = sc_matrix, 
                                         batch_ids = sc_batch, 
                                         cell_type_annotations = sc_celltype_annotations, 
                                         method = method)
###spillover###
simulated_bulk <- simulate_bulk(dataset(annotation = data.frame(ID = colnames(sc_matrix), cell_type = sc_celltype_annotations),
                                        count_matrix = sc_matrix,
                                        name = sc_ds),
                                scaling_factor = "NONE",
                                scenario = "unique", 
                                unique_cell_type="B cell")
signature <- omnideconv::build_model(single_cell_object = as.data.frame(sc_matrix),
                                     cell_type_annotations = sc_celltype_annotations,
                                     bulk_gene_expression = simulated_bulk$pseudo_bulk,
                                     markers = sc_marker,
                                     batch_ids = sc_batch,
                                     verbose=TRUE,
                                     method = method)
deconvolution <- omnideconv::deconvolute(bulk_gene_expression = simulated_bulk$pseudo_bulk,
                                         signature = signature,
                                         single_cell_object = sc_matrix,
                                         batch_ids = sc_batch,
                                         cell_type_annotations = sc_celltype_annotations,
                                         method = method)
