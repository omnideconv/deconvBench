###PREPARE SINGLE CELL DATA AND SAVE INTO RDS###

#'MAYNARD'
#specify resultsdir and author
saveFolder <- "/nfs/data/omnideconv_benchmarking/data/singleCell/maynard"
set <- "maynard_2020_annotated_fine"
#load data
X <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X_tpm.csv", header=FALSE)
sc_tpm <- t(X)
X_counts <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X.csv", header=FALSE)
sc_counts <- t(X_counts)
X_length_scaled <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X_length_scaled.csv", header=FALSE)
sc_length_scaled <- t(X_length_scaled)
obs <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "obs.csv"))
var <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "var.csv"))
colnames(sc_tpm) <- obs$Run
rownames(sc_tpm) <- var$symbol
colnames(sc_counts) <- obs$Run
rownames(sc_counts) <- var$symbol
colnames(sc_length_scaled) <- obs$Run
rownames(sc_length_scaled) <- var$symbol

batch <- obs$patient
marker <- var$symbol

#save data
saveRDS(sc_tpm, file.path(saveFolder, "matrix_tpm.rds"))
saveRDS(sc_counts, file.path(saveFolder, "matrix_counts.rds"))
saveRDS(sc_length_scaled, file.path(saveFolder, "matrix_length_scaled.rds"))
saveRDS(batch, file.path(saveFolder, "batch.rds"))
saveRDS(marker, file.path(saveFolder, "marker.rds"))
saveRDS(obs$cell_type, file.path(saveFolder, "celltype_annotations.rds"))


#'''LAMBRECHTS'''
#specify resultsdir and author
saveFolder <- "/nfs/data/omnideconv_benchmarking/data/singleCell/lambrechts"
set <- "lambrechts_2018_annotated_fine"
#load data
X <- read.csv("/nfs/data/omnideconv_benchmarking/lambrechts_2018_annotated_fine/X_norm_counts.csv", header=FALSE)
sc_norm_counts <- t(X)
X_counts <- read.csv("/nfs/data/omnideconv_benchmarking/lambrechts_2018_annotated_fine/X.csv", header=FALSE)
sc_counts <- t(X_counts)
obs <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "obs.csv"))
var <- read.csv(file.path("/nfs/data/omnideconv_benchmarking", set, "var.csv"))
colnames(sc_norm_counts) <- obs$index
rownames(sc_norm_counts) <- var$X
colnames(sc_counts) <- obs$index
rownames(sc_counts) <- var$X

batch <- obs$patient
marker <- var$X

#save data
saveRDS(sc_norm_counts, file.path(saveFolder, "matrix_norm_counts.rds"))
saveRDS(sc_counts, file.path(saveFolder, "matrix_counts.rds"))
saveRDS(batch, file.path(saveFolder, "batch.rds"))
saveRDS(marker, file.path(saveFolder, "marker.rds"))
saveRDS(obs$cell_type, file.path(saveFolder, "celltype_annotations.rds"))



