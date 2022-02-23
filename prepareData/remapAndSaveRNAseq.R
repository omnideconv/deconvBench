saveDir <- "/nfs/data/omnideconv_benchmarking/data/PBMC"

hoek_pbmc_facs <- read.table("/nfs/data/omnideconv_benchmarking/quanTIseq_data/Hoek_PBMC_FACS.txt", header = TRUE, row.names = 1, sep = "\t") #hoek facs
hoek_pbmc_facs <- data.frame(t(hoek_pbmc_facs))
hoek_mapping <- read.table("/nfs/data/omnideconv_benchmarking/ID_map.txt", header = FALSE, sep = "\t") #hoek mapping
hoek_pbmc_tpm <- read.table("/nfs/data/omnideconv_benchmarking/tpm_Signature_v7.txt", header = TRUE, row.names = 1, sep = "\t") #hoek signature tpm
colnames(hoek_pbmc_tpm) <- hoek_mapping[,1][match(hoek_mapping[,2], colnames(hoek_pbmc_tpm))]
hoek_pbmc_counts <- read.table("/nfs/data/omnideconv_benchmarking/raw_count_Signature_v7.txt", header = TRUE, row.names = 1, sep = "\t") #hoek signature tpm
colnames(hoek_pbmc_counts) <- hoek_mapping[,1][match(hoek_mapping[,2], colnames(hoek_pbmc_counts))]
dataName <- "hoek"
saveRDS(hoek_pbmc_facs, file.path(saveDir, dataName, "hoek_pbmc_facs.rds"))
saveRDS(hoek_pbmc_tpm, file.path(saveDir, dataName, "hoek_pbmc_tpm.rds"))
saveRDS(hoek_pbmc_counts, file.path(saveDir, dataName, "hoek_pbmc_counts.rds"))
#saveRDS(hoek_pbmc_facs, hoek_pbmc_counts, hoek_pbmc_tpm, file.path(saveDir, dataName, "hoek_pbmc_allData.rds"))


finotello_pbmc_tpm <- read.table("/nfs/data/omnideconv_benchmarking/tpm_PBMC_RNAseq.txt", header = TRUE, row.names = 1, sep = "\t") #finotello rnaseq tpm
colnames(finotello_pbmc_tpm) <- gsub("_merged.*", "", colnames(finotello_pbmc_tpm))
finotello_pbmc_counts <- read.table("/nfs/data/omnideconv_benchmarking/raw_count_PBMC_RNAseq.txt", header = TRUE, row.names = 1, sep = "\t") #finotello rnaseq tpm
colnames(finotello_pbmc_counts) <- gsub("_merged.*", "", colnames(finotello_pbmc_counts))
finotello_pbmc_facs <- read.table("/nfs/data/omnideconv_benchmarking/Decon_FACS_20170523.txt", header = TRUE, row.names = 1, sep = "\t")/100 #finotello facs but make it fractions
colnames(finotello_pbmc_facs)<-gsub("D[0]*", "pbmc_", colnames(finotello_pbmc_facs))
rownames(finotello_pbmc_facs) <- gsub("CD8.T.cells", "T.cells.CD8", gsub("CD4.T.cells", "T.cells.CD4", rownames(finotello_pbmc_facs)))
dataName <- "finotello"
saveRDS(finotello_pbmc_facs, file.path(saveDir, dataName, "finotello_pbmc_facs.rds"))
saveRDS(finotello_pbmc_tpm, file.path(saveDir, dataName, "finotello_pbmc_tpm.rds"))
saveRDS(finotello_pbmc_counts, file.path(saveDir, dataName, "finotello_pbmc_counts.rds"))
#saveRDS(finotello_pbmc_facs, finotello_pbmc_counts, finotello_pbmc_tpm, file.path(saveDir, dataName, "finotello_pbmc_allData.rds"))

#Vanderbilt is in one rdata object - counts, gtruth, tpm
load("/nfs/data/omnideconv_benchmarking/vanderbilt/Vanderbilt_Lungcancer.rdata")
vanderbiltLungcancer_pbmc_tpm <- tpm #vanderbilt lungcancer rnaseq tpm
vanderbiltLungcancer_pbmc_counts <- counts #vanderbilt lungcancer rnaseq tpm
vanderbiltLungcancer_pbmc_facs <- gtruth #vanderbilt lungcancer rnaseq facs
dataName <- "vanderbiltLungcancer"
saveRDS(vanderbiltLungcancer_pbmc_facs, file.path(saveDir, dataName, "vanderbiltLungcancer_pbmc_facs.rds"))
saveRDS(vanderbiltLungcancer_pbmc_tpm, file.path(saveDir, dataName, "vanderbiltLungcancer_pbmc_tpm.rds"))
saveRDS(vanderbiltLungcancer_pbmc_counts, file.path(saveDir, dataName, "vanderbiltLungcancer_pbmc_counts.rds"))

load("/nfs/data/omnideconv_benchmarking/vanderbilt/Vanderbilt_Melanoma.rdata")
vanderbiltMelanoma_pbmc_tpm <- tpm #vanderbilt melanoma rnaseq tpm
vanderbiltMelanoma_pbmc_counts <- counts #vanderbilt melanoma rnaseq tpm
vanderbiltMelanoma_pbmc_facs <- gtruth #vanderbilt melanoma rnaseq facs
dataName <- "vanderbiltMelanoma"
saveRDS(vanderbiltMelanoma_pbmc_facs, file.path(saveDir, dataName, "vanderbiltMelanoma_pbmc_facs.rds"))
saveRDS(vanderbiltMelanoma_pbmc_tpm, file.path(saveDir, dataName, "vanderbiltMelanoma_pbmc_tpm.rds"))
saveRDS(vanderbiltMelanoma_pbmc_counts, file.path(saveDir, dataName, "vanderbiltMelanoma_pbmc_counts.rds"))
