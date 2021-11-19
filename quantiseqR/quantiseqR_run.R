### Script that runs immundeconv::quantiseq and quantiseqr with all possible parameters on one RNA-seq data set and stores all results in a RData object ###
library("MASS")

load("/nfs/data/omnideconv_benchmarking/data/PBMC/hoek/hoek_pbmc_tpm.RData") #replace path to RData Object

expression = hoek_pbmc_tpm # change name of matrix in RData object above
dataSetName <- "hoek" # change name of dataset (important for naming results)


for(method in c("lsei", "hampel", "huber", "bisquare")){
	for(ar in c(TRUE, FALSE)){
		for(tu in c(TRUE, FALSE)){
			for(scale in c(TRUE, FALSE)){
				name <- paste(substr(ar, 1, 1), substr(tu, 1, 1), substr(scale, 1, 1), "_", method, sep="")
				quantiseq <- immunedeconv::deconvolute_quantiseq(gene_expression_matrix = expression, arrays = ar, tumor = tu, scale_mrna = scale, method=method) %>% t %>% as.data.frame %>% rownames_to_column(var="Sample") #transformation and extra column so that results are comparable
				assign(paste("quantiseq_", name, "_", dataSetName, sep=""), quantiseq)
				quantiseqR <- quantiseqr::run_quantiseq(expression_data = expression, is_arraydata = ar, is_tumordata = tu, scale_mRNA = scale, method = method)
				assign(paste("quantiseqr_", name, "_", dataSetName, sep=""), quantiseqR)
      }
    }
  }
}

rm(ar)
rm(dataSetName)
rm(method)
rm(name)
rm(scale)
rm(tu)
rm(quantiseq)
rm(quantiseqR)
rm(expression)
rm(hoek_pbmc_tpm)
save.image("/nfs/data/omnideconv_benchmarking/quantiseqComparison/results_hoek.RData") #change path to RData object in which all results should be stored
