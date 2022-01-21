library(ggplot2)

"Usage: 
  umapCelltypes.R <sc_path> <sc_datasetname> <remapping_sheet>
Options:
<sc_path> path to sc dataset
<sc_datasetname> name of sc dataset
<remapping_sheet> excel mapping sheet with remapping" -> doc

print(doc)
args <- docopt::docopt(doc)

sc_path <- args$sc_path
sc_ds <- args$sc_datasetname
sc_matrix <- readRDS(file.path(sc_path, sc_ds, "matrix.rds"))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
remapping_sheet <- args$remapping_sheet

source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)
u <- umap::umap(t(as.matrix(sc_matrix)))
df <- data.frame(u$layout)
df$celltype <- sc_celltype_annotations
p <- ggplot(data = df, aes(X1, X2, color=celltype)) + geom_point() + labs(title=paste("UMAP", sc_ds))
p
ggsave(paste("umap_singlecell_", sc_ds, ".jpg", sep=""), p)

##maynard: looks almost good
##lambrechts: many NAs have to work on mapping again!