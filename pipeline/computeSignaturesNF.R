options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

path <- args[2]
scData <- paste(path, args[4], ".Rdata") #figure out best way to flexibly load data
#currently thought would be RData Object with sc, batch, celltype_annotations
rnaseq_ds <- args[6]
method <- args[8]
outpath <- args[10]

###HARDCODING###
load("/nfs/data/omnideconv_benchmarking/data/PBMC/hoek/hoek_pbmc_tpm.RData")
sc <- t(read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X_tpm.csv", header = FALSE))
obs <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/obs.csv")
var <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/var.csv")
colnames(sc) <- obs$Run
rownames(sc) <- var$symbol
batch <- obs$patient
marker <- var$symbol
################


signature <- omnideconv::build_model(single_cell_object = as.data.frame(sc), cell_type_annotations = obs$cell_type, 
				                                          bulk_gene_expression = hoek_pbmc_tpm, batch_ids = batch, method = "bisque")
# assign(paste("signature_", method, "_", scData, "_", rnaseq_ds, sep=""), signature) HOW DO I SAVE THE RDATA OBJECT THEN??? nvm just save as csv
write.csv(signature, paste(outpath, "/signature_", method, "_", scData, "_", rnaseq_ds, ".csv", sep=""))


