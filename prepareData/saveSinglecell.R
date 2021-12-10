###PREPARE SINGLE CELL DATA AND SAVE INTO RDS###

#'MAYNARD'
#specify resultsdir and author
saveFolder <- "/nfs/data/omnideconv_benchmarking/data/singleCell/maynard"
set <- "maynard_2020_annotated_fine"
#load data
X <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X_tpm.csv", header=FALSE)
sc <- t(X)
obs <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "obs.csv"))
var <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "var.csv"))
colnames(sc) <- obs$Run
rownames(sc) <- var$symbol

batch <- obs$patient
marker <- var$symbol

#save data
saveRDS(sc, file.path(saveFolder, "matrix.rds"))
saveRDS(batch, file.path(saveFolder, "batch.rds"))
saveRDS(marker, file.path(saveFolder, "marker.rds"))
saveRDS(obs$cell_type, file.path(saveFolder, "celltype_annotations.rds"))


#'''LAMBRECHTS'''
#specify resultsdir and author
saveFolder <- "/nfs/data/omnideconv_benchmarking/data/singleCell/lambrechts"
set <- "lambrechts_2018_annotated_fine"
#load data
X <- read.csv("/nfs/data/omnideconv_benchmarking/lambrechts_2018_annotated_fine/X.csv", header = FALSE)
sc <- t(X)
obs <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "obs.csv"))
var <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "var.csv"))
colnames(sc) <- obs$index
rownames(sc) <- var$X


batch <- obs$patient
marker <- var$X

#save data
saveRDS(sc, file.path(saveFolder, "matrix.rds"))
saveRDS(batch, file.path(saveFolder, "batch.rds"))
saveRDS(marker, file.path(saveFolder, "marker.rds"))
saveRDS(obs$cell_type, file.path(saveFolder, "celltype_annotations.rds"))



