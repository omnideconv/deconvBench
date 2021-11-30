options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

path <- args[2]
sc_ds <- args[4]
#load(paste(path, "/", args[4], ".Rdata", sep="")) #figure out best way to flexibly load data
#currently thought would be RData Object with sc, batch, celltype_annotations
rnaseq_ds <- args[6]
load(paste(path, "/", rnaseq_ds, "/", rnaseq_ds, "_pbmc_tpm.RData", sep=""))
method <- args[8]

signature <- readRDS(paste("signature_", method, "_", sc_ds, "_", rnaseq_ds, 
                           ".rds", sep=""))

###HARDCODING###
#load("/nfs/data/omnideconv_benchmarking/data/PBMC/hoek/hoek_pbmc_tpm.RData")
sc <- t(read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X_tpm.csv", header = FALSE))
obs <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/obs.csv")
var <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/var.csv")
colnames(sc) <- obs$Run
rownames(sc) <- var$symbol
batch <- obs$patient
marker <- var$symbol
################

deconvolution <- omnideconv::deconvolute(bulk_gene_expression = get(paste(rnaseq_ds, "_pbmc_tpm", sep="")), 
                        signature = signature, 
                        single_cell_object = sc, 
                        batch_ids = batch, 
                        cell_type_annotations = obs$cell_type, 
                        method = method)

saveRDS(deconvolution, 
        paste("deconvolution_", method, "_", scData, "_", rnaseq_ds, ".rds", sep=""), 
        compress = FALSE)
