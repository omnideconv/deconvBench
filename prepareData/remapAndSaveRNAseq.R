saveDir <- "/nfs/data/omnideconv_benchmarking/data/PBMC"

hoek_pbmc_facs <- read.table("/nfs/data/omnideconv_benchmarking/quanTIseq_data/Hoek_PBMC_FACS.txt", header = TRUE, row.names = 1, sep = "\t") #hoek facs
hoek_pbmc_facs <- data.frame(t(hoek_pbmc_facs))
hoek_mapping <- read.table("/nfs/data/omnideconv_benchmarking/ID_map.txt", header = FALSE, sep = "\t") #hoek mapping
hoek_pbmc_tpm <- read.table("/nfs/data/omnideconv_benchmarking/tpm_Signature_v7.txt", header = TRUE, row.names = 1, sep = "\t") #hoek signature tpm
colnames(hoek_pbmc_tpm) <- hoek_mapping[,1][match(hoek_mapping[,2], colnames(hoek_pbmc_tpm))]
hoek_pbmc_counts <- read.table("/nfs/data/omnideconv_benchmarking/raw_count_Signature_v7.txt", header = TRUE, row.names = 1, sep = "\t") #hoek signature tpm
colnames(hoek_pbmc_counts) <- hoek_mapping[,1][match(hoek_mapping[,2], colnames(hoek_pbmc_counts))]
dataName <- "hoek"
save(hoek_pbmc_facs, file=paste(saveDir, dataName, "hoek_pbmc_facs.RData", sep="/"))
save(hoek_pbmc_tpm, file=paste(saveDir, dataName, "hoek_pbmc_tpm.RData", sep="/"))
save(hoek_pbmc_counts, file=paste(saveDir, dataName, "hoek_pbmc_counts.RData", sep="/"))
save(hoek_pbmc_facs, hoek_pbmc_counts, hoek_pbmc_tpm, file=paste(saveDir, dataName, "hoek_pbmc_allData.RData", sep="/"))


finotello_pbmc_tpm <- read.table("/nfs/data/omnideconv_benchmarking/tpm_PBMC_RNAseq.txt", header = TRUE, row.names = 1, sep = "\t") #finotello rnaseq tpm
colnames(finotello_pbmc_tpm) <- gsub("_merged.*", "", colnames(finotello_pbmc_tpm))
finotello_pbmc_counts <- read.table("/nfs/data/omnideconv_benchmarking/raw_count_PBMC_RNAseq.txt", header = TRUE, row.names = 1, sep = "\t") #finotello rnaseq tpm
colnames(finotello_pbmc_counts) <- gsub("_merged.*", "", colnames(finotello_pbmc_counts))
finotello_pbmc_facs <- read.table("/nfs/data/omnideconv_benchmarking/Decon_FACS_20170523.txt", header = TRUE, row.names = 1, sep = "\t")/100 #finotello facs but make it fractions
colnames(finotello_pbmc_facs)<-gsub("D[0]*", "pbmc_", colnames(finotello_pbmc_facs))
rownames(finotello_pbmc_facs) <- gsub("CD8.T.cells", "T.cells.CD8", gsub("CD4.T.cells", "T.cells.CD4", rownames(finotello_pbmc_facs)))
dataName <- "finotello"
save(finotello_pbmc_facs, file=paste(saveDir, dataName, "finotello_pbmc_facs.RData", sep="/"))
save(finotello_pbmc_tpm, file=paste(saveDir, dataName, "finotello_pbmc_tpm.RData", sep="/"))
save(finotello_pbmc_counts, file=paste(saveDir, dataName, "finotello_pbmc_counts.RData", sep="/"))
save(finotello_pbmc_facs, finotello_pbmc_counts, finotello_pbmc_tpm, file=paste(saveDir, dataName, "finotello_pbmc_allData.RData", sep="/"))

