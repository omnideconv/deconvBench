library(omnideconv)

x <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X_tpm.csv", header = FALSE)
obs <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/obs.csv")
var <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/var.csv")
sc <- t(x)
colnames(sc) <- obs$Run
rownames(sc) <- var$symbol
batch <- obs$patient
marker <- var$symbol
load("/nfs/data/omnideconv_benchmarking/data/PBMC/hoek/hoek_pbmc_tpm.RData")
load("/nfs/data/omnideconv_benchmarking/data/PBMC/hoek/hoek_pbmc_counts.RData")

print("loaded everything")

try({
  sig_bisque <- readRDS("/nfs/data/omnideconv_benchmarking/sig_bisque.RData")
  deconv_bisque <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_bisque, single_cell_object = sc, batch_ids = batch, cell_type_annotations = obs$cell_type, method = "bisque")
  save(deconv_bisque, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_bisque.RData")
})
print("bisque done")
try({sig_autoG <- "/nfs/data/omnideconv_benchmarking/sig_autogenes.pickle"
deconv_autoG <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_autoG, method = "autogenes")
save(deconv_autoG, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_autoG.RData")
})
print("autogenes done")
try({
  sig_bseqsc <- readRDS("/nfs/data/omnideconv_benchmarking/sig_bseqsc.RData")
  deconv_bseqsc <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_bseqsc, method = "bseqsc")
save(deconv_bseqsc, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_bseqsc.RData")
})
print("bseqsc done")
try({deconv_cdseq <- omnideconv::deconvolute_cdseq(bulk_gene_expression = hoek_pbmc_tpm, single_cell_object = sc, batch_ids = batch, cell_type_annotations = obs$cell_type)
save(deconv_cdseq, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_cdseq.RData")
})
print("cdseq done")
try({
  sig_cibersortx <- load("/nfs/data/omnideconv_benchmarking/sig_cibersortx.RData")
  deconv_cibersortx <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_cibersortx, method = "cibersortx")
  save(deconv_cibersortx, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_cibersortx.RData")
})
print("cibersortx done")
try({deconv_cpm <- omnideconv::deconvolute_cpm(bulk_gene_expression = hoek_pbmc_tpm, single_cell_object = sc, cell_type_annotations = obs$cell_type)
save(deconv_cpm, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_cpm.RData")
})
print("cpm done")
try({
  sig_dwls <- readRDS("/nfs/data/omnideconv_benchmarking/sig_dwls.RData")
deconv_dwls <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_dwls, method = "dwls")
save(deconv_dwls, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_dwls.RData")
})
print("dwls done")
try({
  load("/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_momf.RData")
deconv_momf <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_momf, single_cell_object = sc, method = "momf")
save(deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_momf, single_cell_object = sc, method = "momf"), file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_momf.RData")
})
print("momf done")
try({
deconv_music <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = NULL, single_cell_object = sc, batch_ids = batch, cell_type_annotations = obs$cell_type, method="music")
save(deconv_music, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_music.RData")
})
print("music done")
try({
  sig_scaden <- load("/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_scaden.RData")
deconv_scaden <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = sig_scaden, method = "scaden")
save(deconv_scaden, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_scaden.RData")
})
print("scaden done")
try({
deconv_scdc <- omnideconv::deconvolute(bulk_gene_expression = hoek_pbmc_tpm, signature = NULL, single_cell_object = sc, batch_ids = batch, cell_type_annotations = obs$cell_type, method="scdc")
save(deconv_scdc, file="/nfs/data/omnideconv_benchmarking/deconv_hoek_gtruth/deconv_scdc.RData")
})
print("scdc done")

