options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

path <- args[2]
sc_ds <- args[4]
sc_path <- paste(path, args[4], ".Rdata", sep="") #figure out best way to flexibly load data
#currently thought would be RData Object with sc, batch, celltype_annotations
rnaseq_ds <- args[6]
load(paste(path, "/", rnaseq_ds, "/", rnaseq_ds, "_pbmc_tpm.RData", sep=""))
method <- args[8]

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
if(method=="cibersortx"){
  omnideconv::set_cibersortx_credentials("k.reinisch@campus.lmu.de", "a05114832330fda42ce0f5596875ee0d")
}

signature <- omnideconv::build_model(single_cell_object = as.data.frame(sc), 
                                     cell_type_annotations = obs$cell_type, 
				                            bulk_gene_expression = get(paste(rnaseq_ds, "_pbmc_tpm", sep="")), 
				                            batch_ids = batch, 
				                            method = method)
#if method has no build_model, what happens?
saveRDS(signature, 
        paste("signature_", method, "_", sc_ds, "_", rnaseq_ds, ".rds", sep=""), 
        compress = FALSE)

