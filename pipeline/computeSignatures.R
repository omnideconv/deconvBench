library("omnideconv")

x <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/X_tpm.csv", header = FALSE)
obs <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/obs.csv")
var <- read.csv("/nfs/data/omnideconv_benchmarking/maynard_2020_annotated_fine/var.csv")
load("/nfs/data/omnideconv_benchmarking/hoek_sample_annotations.RData") #loads mix.mat
sc <- t(x)
colnames(sc) <- obs$Run
rownames(sc) <- var$symbol
batch <- obs$patient
marker <- var$symbol

try({sig_autoG <- omnideconv::build_model(as.data.frame(sc), obs$cell_type, method = "autogenes")
save(sig_autoG, file="/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_autogenes.RData") #only a path
})
try({sig_bisque <- omnideconv::build_model(as.data.frame(sc), obs$cell_type, batch_ids = batch, method = "bisque")
save(sig_bisque, file="/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_bisque.RData")
})
try({
  sig_bseqsc <- omnideconv::build_model(single_cell_object = sc, cell_type_annotations = obs$cell_type, batch_ids = colnames(sc), markers = marker, limit=TRUE, method = "bseqsc")
save(sig_bseqsc, file="/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_bseqsc.RData")
})
#sig_cdseq <- build_model(as.data.frame(sc), obs$cell_type, method = "cdseq") #only in one step
try({set_cibersortx_credentials("k.reinisch@campus.lmu.de", "a05114832330fda42ce0f5596875ee0d")
  sig_cibersortx <- omnideconv::build_model(as.data.frame(sc), obs$cell_type, outputdir="/nfs/data/omnideconv_benchmarking/data/maynard_sig" , method = "cibersortx")
  save(sig_cibersortx, file="/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_cibersortx.RData")
})
#sig_cpm <- build_model(as.data.frame(sc), obs$cell_type, bulk_gene_expression = hoek_star_rsem, method = "cpm") #only in one step
try({
  sig_dwls <- omnideconv::build_model(as.data.frame(sc), obs$cell_type, method = "dwls")
  save(sig_dwls, file="/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_dwls.RData")
})
try({sig_momf <- omnideconv::build_model(as.data.frame(sc), obs$cell_type, bulk_gene_expression = hoek_pbmc_tpm, method = "momf")
save(sig_momf, file="/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_momf.RData")
})
#sig_music <- build_model(as.data.frame(sc), obs$cell_type, batch_ids = batch, method = "music") #only in one step
try({sig_scaden <- omnideconv::build_model(as.data.frame(sc), obs$cell_type, bulk_gene_expression = hoek_pbmc_tpm, method = "scaden")
save(sig_scaden, file="/nfs/data/omnideconv_benchmarking/data/maynard_sig/sig_scaden.RData") #only a path
})
#sig_scdc <- build_model(as.data.frame(sc), obs$cell_type, batch_ids = batch, method = "scdc") #only in one step

