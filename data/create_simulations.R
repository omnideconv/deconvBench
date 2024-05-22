library(tidyverse)
library(SimBu)
library(data.table)
library(parallel)

set.seed(123)

create_simulation <- function(bulk_features, simbu_ds){
  
  res <- parallel::mclapply(1:bulk_features$nsamples, function(i){
    print(paste0('Simulating technical replicates based on real sample ',names(bulk_features$read_counts)[i],' ...'))
    
    custom_fractions <- bulk_features$proportions[,i]
    names(custom_fractions) <- rownames(bulk_features$proportions)
    # only use cell types that are also present in reference
    custom_fractions <- custom_fractions[names(custom_fractions) %in% unique(simbu_ds@colData$cell_type)]
    # scale custom fractions to sum up to 1
    custom_fractions <- custom_fractions/sum(custom_fractions)
    # copy these fractions 10 times for 10 replicates
    custom_fractions_df <- data.frame(matrix(rep(custom_fractions, 10), nrow = 10, byrow = TRUE))
    colnames(custom_fractions_df) <- names(custom_fractions)
    rownames(custom_fractions_df) <- paste0('custom_sample',1:10)
    
    print('Simulating with bias ...')
    simulation <- SimBu::simulate_bulk(simbu_ds, 
                                       scenario = 'custom', 
                                       scaling_factor = 'expressed_genes', 
                                       custom_scenario_data = custom_fractions_df, 
                                       nsamples = 10, 
                                       ncells = 10000, 
                                       total_read_counts = bulk_features$read_counts[i], 
                                       BPPARAM = BiocParallel::MulticoreParam(workers=10), 
                                       run_parallel=F)
    
    print('Simulating without bias ...')
    simulation_nobias <- SimBu::simulate_bulk(simbu_ds, 
                                              scenario = 'custom', 
                                              scaling_factor = 'NONE', 
                                              custom_scenario_data = custom_fractions_df, 
                                              nsamples = 10, 
                                              ncells = 10000, 
                                              total_read_counts = bulk_features$read_counts[i], 
                                              BPPARAM = BiocParallel::MulticoreParam(workers=10), 
                                              run_parallel=F)
    
    return(list('bias'=simulation, 'nobias'=simulation_nobias))
    
  }, mc.cores = bulk_features$nsamples)
  
  simulation_full <- SimBu::merge_simulations(lapply(res, function(x) x$bias))
  simulation_nobias_full <- SimBu::merge_simulations(lapply(res, function(x) x$nobias))
  
  return(list('bias'=simulation_full, 'nobias'=simulation_nobias_full))
}


#### Mouse

chen_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen/chen_counts.rds')
chen_facs <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen/chen_facs.rds')
chen_features <- list('read_counts'=colSums(chen_counts), 'proportions'=chen_facs[,colnames(chen_counts)], 'nsamples'=dim(chen_counts)[2])

peti_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez/hoek_counts.rds')
peti_facs <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez/hoek_facs.rds')
peti_features <- list('read_counts'=colSums(peti_counts), 'proportions'=peti_facs[,colnames(peti_counts)], 'nsamples'=dim(peti_counts)[2])

tm_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/matrix_counts.rds')
tm_norm_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula/matrix_norm_counts.rds')
tm_anno <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula/celltype_annotations.rds')
tm_anno_df <- data.frame('ID'=colnames(tm_counts), 'cell_type'=tm_anno)

tm_ds <- SimBu::dataset(annotation = tm_anno_df, 
                         count_matrix = tm_counts, 
                         tpm_matrix = tm_norm_counts, 
                         name = 'tabula-muris', 
                         filter_genes = T)

## Chen based simulations
chen_simulations <- create_simulation(chen_features, tm_ds)

colnames(chen_simulations$bias$bulk) <- rownames(chen_simulations$bias$cell_fractions)
saveRDS(SummarizedExperiment::assays(chen_simulations$bias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation/chen-simulation_counts.rds')
saveRDS(SummarizedExperiment::assays(chen_simulations$bias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation/chen-simulation_tpm.rds')
saveRDS(t(as.matrix(chen_simulations$bias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation/chen-simulation_facs.rds')

colnames(chen_simulations$nobias$bulk) <- rownames(chen_simulations$nobias$cell_fractions)
saveRDS(SummarizedExperiment::assays(chen_simulations$nobias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation-nobias/chen-simulation-nobias_counts.rds')
saveRDS(SummarizedExperiment::assays(chen_simulations$nobias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation-nobias/chen-simulation-nobias_tpm.rds')
saveRDS(t(as.matrix(chen_simulations$nobias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation-nobias/chen-simulation-nobias_facs.rds')


## Petitprez based simulations
peti_simulations <- create_simulation(peti_features, tm_ds)

colnames(peti_simulations$bias$bulk) <- rownames(peti_simulations$bias$cell_fractions)
saveRDS(SummarizedExperiment::assays(peti_simulations$bias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation/petitprez-simulation_counts.rds')
saveRDS(SummarizedExperiment::assays(peti_simulations$bias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation/petitprez-simulation_tpm.rds')
saveRDS(t(as.matrix(peti_simulations$bias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation/petitprez-simulation_facs.rds')

colnames(peti_simulations$nobias$bulk) <- rownames(peti_simulations$nobias$cell_fractions)
saveRDS(SummarizedExperiment::assays(peti_simulations$nobias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation-nobias/petitprez-simulation-nobias_counts.rds')
saveRDS(SummarizedExperiment::assays(peti_simulations$nobias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation-nobias/petitprez-simulation-nobias_tpm.rds')
saveRDS(t(as.matrix(peti_simulations$nobias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation-nobias/petitprez-simulation-nobias_facs.rds')

#### Human

fino_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello/finotello_counts.rds')
fino_facs <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello/finotello_facs.rds')
fino_features <- list('read_counts'=colSums(fino_counts), 'proportions'=fino_facs[,colnames(fino_counts)], 'nsamples'=dim(fino_counts)[2])

hoek_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek/hoek_counts.rds')
hoek_facs <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek/hoek_facs.rds')
hoek_features <- list('read_counts'=colSums(hoek_counts), 'proportions'=hoek_facs[,colnames(hoek_counts)], 'nsamples'=dim(hoek_counts)[2])

hao_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/matrix_counts_dense.rds')
hao_norm_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/matrix_norm_counts_dense.rds')
hao_anno <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/hao-complete/celltype_annotations.rds')
hao_anno_df <- data.frame('ID'=colnames(hao_counts), 'cell_type'=hao_anno)

hao_ds <- SimBu::dataset(annotation = hao_anno_df, 
                         count_matrix = hao_counts, 
                         tpm_matrix = hao_norm_counts, 
                         name = 'hao', 
                         filter_genes = T)

## finotello based simulations
fino_simulations <- create_simulation(fino_features, hao_ds)

colnames(fino_simulations$bias$bulk) <- rownames(fino_simulations$bias$cell_fractions)
saveRDS(SummarizedExperiment::assays(fino_simulations$bias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation/finotello-simulation_counts.rds')
saveRDS(SummarizedExperiment::assays(fino_simulations$bias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation/finotello-simulation_tpm.rds')
saveRDS(t(as.matrix(fino_simulations$bias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation/finotello-simulation_facs.rds')

colnames(fino_simulations$nobias$bulk) <- rownames(fino_simulations$nobias$cell_fractions)
saveRDS(SummarizedExperiment::assays(fino_simulations$nobias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation-nobias/finotello-simulation-nobias_counts.rds')
saveRDS(SummarizedExperiment::assays(fino_simulations$nobias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation-nobias/finotello-simulation-nobias_tpm.rds')
saveRDS(t(as.matrix(fino_simulations$nobias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation-nobias/finotello-simulation-nobias_facs.rds')


## hoek based simulations

# change annotation in Hao of all T cell subtypes to regular T cells, so that we can sample them to match Hoek FACS
hao_anno_df_tcells <- hao_anno_df %>% 
  mutate(cell_type = recode(cell_type,
                            'T cells CD4 conv'='T cell',
                            'T cells CD8'='T cell',
                            'Tregs'='T cell'))

hao_ds_tcells <- SimBu::dataset(annotation = hao_anno_df_tcells, 
                                count_matrix = hao_counts, 
                                tpm_matrix = hao_norm_counts, 
                                name = 'hao', 
                                filter_genes = T)

hoek_simulations <- create_simulation(hoek_features, hao_ds_tcells)

colnames(hoek_simulations$bias$bulk) <- rownames(hoek_simulations$bias$cell_fractions)
saveRDS(SummarizedExperiment::assays(hoek_simulations$bias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation/hoek-simulation_counts.rds')
saveRDS(SummarizedExperiment::assays(hoek_simulations$bias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation/hoek-simulation_tpm.rds')
saveRDS(t(as.matrix(hoek_simulations$bias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation/hoek-simulation_facs.rds')

colnames(hoek_simulations$nobias$bulk) <- rownames(hoek_simulations$nobias$cell_fractions)
saveRDS(SummarizedExperiment::assays(hoek_simulations$nobias$bulk)[["bulk_counts"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation-nobias/hoek-simulation-nobias_counts.rds')
saveRDS(SummarizedExperiment::assays(hoek_simulations$nobias$bulk)[["bulk_tpm"]], 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation-nobias/hoek-simulation-nobias_tpm.rds')
saveRDS(t(as.matrix(hoek_simulations$nobias$cell_fractions)), 
        '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation-nobias/hoek-simulation-nobias_facs.rds')

