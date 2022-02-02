#!/usr/bin/Rscript

"Usage: 
  deconvoluteSimBulkNF.R <sc_path> <sc_datasetname> <bulk> <deconv_method> <remapping_sheet>
Options:
<sc_path> path to sc dataset
<sc_datasetname> name of sc dataset
<bulk> simulated bulk
<deconv_method>  deconv method
<remapping_sheet> excel mapping sheet with remapping" -> doc

args <- docopt::docopt(doc)

sc_path <- args$sc_path
sc_ds <- args$sc_datasetname
sc_matrix <- readRDS(file.path(sc_path, sc_ds, "matrix.rds"))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))

simulated_bulk <- readRDS(args$bulk)
bulk_scenario <- gsub("simulatedBulk_", "", gsub(".rds", "", basename(args$bulk)))

method <- args$deconv_method
remapping_sheet <- args$remapping_sheet

tissue <- args$tissue
organism <- args$organism

source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)

if(method=="cibersortx"){
  omnideconv::set_cibersortx_credentials("k.reinisch@campus.lmu.de", 
                                         "a05114832330fda42ce0f5596875ee0d")
}

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
saveRDS(deconvolution, 
        paste("deconvolution_", bulk_scenario, "_", method, ".rds", sep=""), 
        compress = FALSE)


if(grepl("spillover", bulk_scenario)){
  names <- c("B cell")
  pdf(paste("empty.pdf", sep=""))
  
  dev.off()
} else {
  names <- c("Monocyte conventional", "Macrophage", "Endothelial cell",
             "Mast cell", "B cell", "T cell CD8+", "T cell CD4+", "Epithelial cell",
             "other cell", "NK cell", "Monocyte non-conventional",
             "Myeloid dendritic cell", "Cancer associated fibroblast", "Club cell",
             "B cell plasma", "T cell regulatory (Tregs)", "Plasmacytoid dendritic cell",
             "Langerhans cell", "Goblet cell")
  colnames(simulated_bulk$cell_fractions) <- names
  corMatrix <- cor(deconvolution, simulated_bulk$cell_fractions)
  pdf(paste("correlationMatrix_", bulk_scenario, "_", method, ".pdf", sep=""))
  corrplot::corrplot(corMatrix, 
                   is.corr = TRUE, tl.col = "black", order = 'alphabet', 
                   title = paste(bulk_scenario, "-", method))
  dev.off()
}