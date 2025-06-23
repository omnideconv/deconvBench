library(SimBu)
library(dplyr)
library(Matrix)
library(SummarizedExperiment)
set.seed(123)

#### Mouse ####
tabula_muris_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/matrix_counts.rds')
tabula_muris_norm_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/matrix_norm_counts.rds')
tabula_muris_annotation <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tabula-muris/celltype_annotations.rds')

tabula_muris_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(tabula_muris_counts), cell_type = tabula_muris_annotation), 
  count_matrix = tabula_muris_counts,
  tpm_matrix = tabula_muris_norm_counts,
  name = 'tabula-muris', 
  filter_genes = T)

chen_bulk_celltypes <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen/chen_facs.rds')
chen_bulk_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen/chen_counts.rds')
petitprez_bulk_celltypes <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez/petitprez_facs.rds')
petitprez_bulk_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez/petitprez_counts.rds')

sc_celltypes <- unique(tabula_muris_annotation)

chen_matching_celltypes <- chen_bulk_celltypes[sc_celltypes,] %>% 
  subset(!grepl('NA', row.names(.))) %>% 
  select_if(~ !any(is.na(.))) %>%
  scale(., center=F, scale=colSums(.))
petitprez_matching_celltypes <- petitprez_bulk_celltypes[sc_celltypes,] %>% 
  subset(!grepl('NA', row.names(.))) %>% 
  select_if(~ !any(is.na(.))) %>%
  scale(., center=F, scale=colSums(.))

##### create simulations #####
chen_seed <- 1234
sim_list_chen <- lapply(1:ncol(chen_matching_celltypes), function(i){
  sample_id <- colnames(chen_matching_celltypes)[i]
  total_mapped_reads <- colSums(chen_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = tabula_muris_ds, 
                       scenario = 'custom', 
                       scaling_factor = 'expressed_genes', 
                       nsamples = 10, 
                       ncells = 10000, 
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, chen_matching_celltypes[,i])), check.names = F), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10), 
                       run_parallel = T, 
                       seed = chen_seed)
})
sim_chen <- SimBu::merge_simulations(sim_list_chen)
colnames(sim_chen$bulk) <- paste0(colnames(sim_chen$bulk), "_", seq_len(length(colnames(sim_chen$bulk))))
saveRDS(assays(sim_chen$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation/chen-simulation_counts.rds')
saveRDS(assays(sim_chen$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation/chen-simulation_tpm.rds')
saveRDS(t(sim_chen$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation/chen-simulation_facs.rds')

sim_list_chen_noBias <- lapply(1:ncol(chen_matching_celltypes), function(i){
  sample_id <- colnames(chen_matching_celltypes)[i]
  total_mapped_reads <- colSums(chen_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = tabula_muris_ds, 
                       scenario = 'custom', 
                       scaling_factor = 'NONE', 
                       nsamples = 10, 
                       ncells = 10000, 
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, chen_matching_celltypes[,i])), check.names = F), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10), 
                       run_parallel = T, 
                       seed = chen_seed)
})
sim_chen_noBias <- SimBu::merge_simulations(sim_list_chen_noBias)
colnames(sim_chen_noBias$bulk) <- paste0(colnames(sim_chen_noBias$bulk), "_", seq_len(length(colnames(sim_chen_noBias$bulk))))
saveRDS(assays(sim_chen_noBias$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation-nobias/chen-simulation-nobias_counts.rds')
saveRDS(assays(sim_chen_noBias$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation-nobias/chen-simulation-nobias_tpm.rds')
saveRDS(t(sim_chen_noBias$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/chen-simulation-nobias/chen-simulation-nobias_facs.rds')


petiprez_seed <- 2345
sim_list_petitprez <- lapply(1:ncol(petitprez_matching_celltypes), function(i){
  sample_id <- colnames(petitprez_matching_celltypes)[i]
  total_mapped_reads <- colSums(petitprez_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = tabula_muris_ds, 
                       scenario = 'custom', 
                       scaling_factor = 'expressed_genes', 
                       nsamples = 10, 
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, petitprez_matching_celltypes[,i])), check.names = F), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10), 
                       run_parallel = T, 
                       seed = petiprez_seed)
})
sim_petitprez <- SimBu::merge_simulations(sim_list_petitprez)
colnames(sim_petitprez$bulk) <- paste0(colnames(sim_petitprez$bulk), "_", seq_len(length(colnames(sim_petitprez$bulk))))
saveRDS(assays(sim_petitprez$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation/petitprez-simulation_counts.rds')
saveRDS(assays(sim_petitprez$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation/petitprez-simulation_tpm.rds')
saveRDS(t(sim_petitprez$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation/petitprez-simulation_facs.rds')

sim_list_petitprez_noBias <- lapply(1:ncol(petitprez_matching_celltypes), function(i){
  sample_id <- colnames(petitprez_matching_celltypes)[i]
  total_mapped_reads <- colSums(petitprez_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = tabula_muris_ds, 
                       scenario = 'custom', 
                       scaling_factor = 'NONE', 
                       nsamples = 10, 
                       ncells = 10000, 
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, petitprez_matching_celltypes[,i])), check.names = F), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10), 
                       run_parallel = T, 
                       seed = petiprez_seed)
})
sim_petitprez_noBias <- SimBu::merge_simulations(sim_list_petitprez_noBias)
colnames(sim_petitprez_noBias$bulk) <- paste0(colnames(sim_petitprez_noBias$bulk), "_", seq_len(length(colnames(sim_petitprez_noBias$bulk))))
saveRDS(assays(sim_petitprez_noBias$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation-nobias/petitprez-simulation-nobias_counts.rds')
saveRDS(assays(sim_petitprez_noBias$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation-nobias/petitprez-simulation-nobias_tpm.rds')
saveRDS(t(sim_petitprez_noBias$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/Mouse/petitprez-simulation-nobias/petitprez-simulation-nobias_facs.rds')



#### Human PBMC ####

hao_ds <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCellXGeneCensus/HaoCleaned/simbu_ds.rds')

# create additional levels of annotation for Hoek and Altman
hao_ds_hoek <- hao_ds
hao_ds_hoek@colData$cell_type <- recode(hao_ds_hoek@colData$cell_type, 
                                        'T cells CD8' = 'T cell',
                                        'T cells CD4 conv' = 'T cell',
                                        'Tregs' = 'T cell')
hao_ds_altman <- hao_ds
hao_ds_altman@colData$cell_type <- recode(hao_ds_altman@colData$cell_type, 
                                        'T cells CD8' = 'Lymphocytes',
                                        'T cells CD4 conv' = 'Lymphocytes',
                                        'Tregs' = 'Lymphocytes',
                                        'B cells' = 'Lymphocytes',
                                        'Plasma cells' = 'Lymphocytes',
                                        'NK cells' = 'Lymphocytes',
                                        'ILC' = 'Lymphocytes',
                                        'mDC' = 'Monocytes',
                                        'pDC' = 'Monocytes')


finotello_bulk_celltypes <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello/finotello_facs.rds')
finotello_bulk_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello/finotello_counts.rds')
hoek_bulk_celltypes <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek/hoek_facs.rds')
hoek_bulk_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek/hoek_counts.rds')
morandini_bulk_celltypes <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini/morandini_facs.rds')
morandini_bulk_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini/morandini_counts.rds')
altman_bulk_celltypes <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman/altman_facs.rds')
altman_bulk_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman/altman_counts.rds')

finotello_bulk_celltypes <- finotello_bulk_celltypes[,colnames(finotello_bulk_counts)]
morandini_bulk_counts <- morandini_bulk_counts[,colnames(morandini_bulk_celltypes)]

hoek_matching_celltypes <- hoek_bulk_celltypes %>% 
  subset(row.names(.) %in% unique(hao_ds_hoek@colData$cell_type)) %>%  
  select_if(~ !any(is.na(.))) %>%
  scale(., center=F, scale=colSums(.))

finotello_matching_celltypes <- finotello_bulk_celltypes %>% 
  subset(row.names(.) %in% unique(hao_ds@colData$cell_type)) %>% 
  select_if(~ !any(is.na(.))) %>%
  scale(., center=F, scale=colSums(.))

rownames(morandini_bulk_celltypes)[which(rownames(morandini_bulk_celltypes) == 'T cells CD4')] <- 'T cells CD4 conv'
morandini_matching_celltypes <- morandini_bulk_celltypes %>% 
  subset(row.names(.) %in% unique(hao_ds@colData$cell_type)) %>% 
  select_if(~ !any(is.na(.))) %>%
  scale(., center=F, scale=colSums(.))

altman_matching_celltypes <- altman_bulk_celltypes %>% 
  subset(row.names(.) %in% unique(hao_ds_altman@colData$cell_type)) %>% 
  select_if(~ !any(is.na(.))) %>%
  scale(., center=F, scale=colSums(.))



##### create simulations #####

###### HaoCleaned #####
sc_ds <- 'HaoCleaned'
ncells <- 10000
nsamples <- 50
scenario <- 'random'
bias <- 'expressed_genes'
data_dir_bulk <- '/nfs/data/omnideconv_benchmarking_clean/data/simulations/pseudo_bulk_mrnabias/'

pseudobulk_name <- paste0(sc_ds, "-ncells", ncells, "-nsamples", nsamples, "-", scenario, "-", bias, "-simulation")
output_dir <- paste0(data_dir_bulk, '/', pseudobulk_name)

set.seed(1234)
simulated_bulk <- SimBu::simulate_bulk(
  data =  hao_ds,
  scenario = scenario,
  scaling_factor = bias,
  nsamples = nsamples,
  ncells = ncells,
  total_read_counts = 100000000,  
  BPPARAM = BiocParallel::MulticoreParam(workers = 20),
  run_parallel = TRUE,  
  seed = 1234
)

saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_counts"]], paste0(output_dir,'/', pseudobulk_name, '_counts.rds'))
saveRDS(SummarizedExperiment::assays(simulated_bulk$bulk)[["bulk_tpm"]], paste0(output_dir,'/', pseudobulk_name, '_tpm.rds'))
saveRDS(t(simulated_bulk$cell_fractions), paste0(output_dir,'/', pseudobulk_name, '_facs.rds'))

###### hoek ######
hoek_seed <- 3456
sim_list_hoek <- lapply(1:ncol(hoek_matching_celltypes), function(i){
  sample_id <- colnames(hoek_matching_celltypes)[i]
  total_mapped_reads <- colSums(hoek_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds_hoek,
                       scenario = 'custom',
                       scaling_factor = 'expressed_genes',
                       nsamples = 10,
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, hoek_matching_celltypes[,i])), check.names = F),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10),
                       run_parallel = T, 
                       seed = hoek_seed)
})
sim_hoek <- SimBu::merge_simulations(sim_list_hoek)
colnames(sim_hoek$bulk) <- paste0(colnames(sim_hoek$bulk), "_", seq_len(length(colnames(sim_hoek$bulk))))
saveRDS(assays(sim_hoek$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation/hoek-simulation_counts.rds')
saveRDS(assays(sim_hoek$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation/hoek-simulation_tpm.rds')
saveRDS(t(sim_hoek$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation/hoek-simulation_facs.rds')


sim_list_hoek_noBias <- lapply(1:ncol(hoek_matching_celltypes), function(i){
  sample_id <- colnames(hoek_matching_celltypes)[i]
  total_mapped_reads <- colSums(hoek_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds_hoek,
                       scenario = 'custom',
                       scaling_factor = 'NONE',
                       nsamples = 10,
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, hoek_matching_celltypes[,i])), check.names = F),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10),
                       run_parallel = T, 
                       seed = hoek_seed)
})
sim_hoek_noBias <- SimBu::merge_simulations(sim_list_hoek_noBias)
colnames(sim_hoek_noBias$bulk) <- paste0(colnames(sim_hoek_noBias$bulk), "_", seq_len(length(colnames(sim_hoek_noBias$bulk))))
saveRDS(assays(sim_hoek_noBias$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation-nobias/hoek-simulation-nobias_counts.rds')
saveRDS(assays(sim_hoek_noBias$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation-nobias/hoek-simulation-nobias_tpm.rds')
saveRDS(t(sim_hoek_noBias$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/hoek-simulation-nobias/hoek-simulation-nobias_facs.rds')


###### finotello ######
fino_seed <- 4567
sim_list_finotello <- lapply(1:ncol(finotello_matching_celltypes), function(i){
  sample_id <- colnames(finotello_matching_celltypes)[i]
  total_mapped_reads <- colSums(finotello_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds,
                       scenario = 'custom',
                       scaling_factor = 'expressed_genes',
                       nsamples = 10,
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, finotello_matching_celltypes[,i])), check.names = F),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10),
                       run_parallel = T, 
                       seed = fino_seed)
})
sim_finotello <- SimBu::merge_simulations(sim_list_finotello)
colnames(sim_finotello$bulk) <- paste0(colnames(sim_finotello$bulk), "_", seq_len(length(colnames(sim_finotello$bulk))))
saveRDS(assays(sim_finotello$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation/finotello-simulation_counts.rds')
saveRDS(assays(sim_finotello$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation/finotello-simulation_tpm.rds')
saveRDS(t(sim_finotello$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation/finotello-simulation_facs.rds')


sim_list_finotello_noBias <- lapply(1:ncol(finotello_matching_celltypes), function(i){
  sample_id <- colnames(finotello_matching_celltypes)[i]
  total_mapped_reads <- colSums(finotello_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds, 
                       scenario = 'custom', 
                       scaling_factor = 'NONE', 
                       nsamples = 10, 
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(10, finotello_matching_celltypes[,i])), check.names = F), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 10), 
                       run_parallel = T, 
                       seed = fino_seed)
})
sim_finotello_noBias <- SimBu::merge_simulations(sim_list_finotello_noBias)
colnames(sim_finotello_noBias$bulk) <- paste0(colnames(sim_finotello_noBias$bulk), "_", seq_len(length(colnames(sim_finotello_noBias$bulk))))
saveRDS(assays(sim_finotello_noBias$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation-nobias/finotello-simulation-nobias_counts.rds')
saveRDS(assays(sim_finotello_noBias$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation-nobias/finotello-simulation-nobias_tpm.rds')
saveRDS(t(sim_finotello_noBias$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/finotello-simulation-nobias/finotello-simulation-nobias_facs.rds')

###### morandini ######

mora_seed <- 45678
sim_list_morandini <- lapply(1:ncol(morandini_matching_celltypes), function(i){
  sample_id <- colnames(morandini_matching_celltypes)[i]
  total_mapped_reads <- colSums(morandini_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds,
                       scenario = 'custom',
                       scaling_factor = 'expressed_genes',
                       nsamples = 5,
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(5, morandini_matching_celltypes[,i])), check.names = F),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 5),
                       run_parallel = T, 
                       seed = mora_seed)
})
sim_morandini <- SimBu::merge_simulations(sim_list_morandini)
colnames(sim_morandini$bulk) <- paste0(colnames(sim_morandini$bulk), "_", seq_len(length(colnames(sim_morandini$bulk))))
saveRDS(assays(sim_morandini$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini-simulation/morandini-simulation_counts.rds')
saveRDS(assays(sim_morandini$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini-simulation/morandini-simulation_tpm.rds')
saveRDS(t(sim_morandini$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini-simulation/morandini-simulation_facs.rds')

sim_list_morandini_noBias <- lapply(1:ncol(morandini_matching_celltypes), function(i){
  sample_id <- colnames(morandini_matching_celltypes)[i]
  total_mapped_reads <- colSums(morandini_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds, 
                       scenario = 'custom', 
                       scaling_factor = 'NONE', 
                       nsamples = 5, 
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(5, morandini_matching_celltypes[,i])), check.names = F), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 5), 
                       run_parallel = T, 
                       seed = mora_seed)
})
sim_morandini_noBias <- SimBu::merge_simulations(sim_list_morandini_noBias)
colnames(sim_morandini_noBias$bulk) <- paste0(colnames(sim_morandini_noBias$bulk), "_", seq_len(length(colnames(sim_morandini_noBias$bulk))))
saveRDS(assays(sim_morandini_noBias$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini-simulation-nobias/morandini-simulation-nobias_counts.rds')
saveRDS(assays(sim_morandini_noBias$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini-simulation-nobias/morandini-simulation-nobias_tpm.rds')
saveRDS(t(sim_morandini_noBias$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/morandini-simulation-nobias/morandini-simulation-nobias_facs.rds')


###### altman ######

altm_seed <- 4789
sim_list_altman <- parallel::mclapply(1:ncol(altman_matching_celltypes), function(i){
  sample_id <- colnames(altman_matching_celltypes)[i]
  total_mapped_reads <- colSums(altman_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds_altman,
                       scenario = 'custom',
                       scaling_factor = 'expressed_genes',
                       nsamples = 5,
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(5, altman_matching_celltypes[,i])), check.names = F),
                       BPPARAM = BiocParallel::MulticoreParam(workers = 5),
                       run_parallel = T, 
                       seed = altm_seed)
}, mc.cores=20)
sim_altman <- SimBu::merge_simulations(sim_list_altman)
colnames(sim_altman$bulk) <- paste0(colnames(sim_altman$bulk), "_", seq_len(length(colnames(sim_altman$bulk))))
saveRDS(assays(sim_altman$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman-simulation/altman-simulation_counts.rds')
saveRDS(assays(sim_altman$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman-simulation/altman-simulation_tpm.rds')
saveRDS(t(sim_altman$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman-simulation/altman-simulation_facs.rds')

sim_list_altman_noBias <- parallel::mclapply(1:ncol(altman_matching_celltypes), function(i){
  sample_id <- colnames(altman_matching_celltypes)[i]
  total_mapped_reads <- colSums(altman_bulk_counts)[sample_id]
  SimBu::simulate_bulk(data = hao_ds_altman, 
                       scenario = 'custom', 
                       scaling_factor = 'NONE', 
                       nsamples = 5, 
                       ncells = 10000,
                       total_read_counts = total_mapped_reads,
                       custom_scenario_data = data.frame(t(replicate(5, altman_matching_celltypes[,i])), check.names = F), 
                       BPPARAM = BiocParallel::MulticoreParam(workers = 5), 
                       run_parallel = T, 
                       seed = altm_seed)
}, mc.cores=20)
sim_altman_noBias <- SimBu::merge_simulations(sim_list_altman_noBias)
colnames(sim_altman_noBias$bulk) <- paste0(colnames(sim_altman_noBias$bulk), "_", seq_len(length(colnames(sim_altman_noBias$bulk))))
saveRDS(assays(sim_altman_noBias$bulk)[['bulk_counts']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman-simulation-nobias/altman-simulation-nobias_counts.rds')
saveRDS(assays(sim_altman_noBias$bulk)[['bulk_tpm']], '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman-simulation-nobias/altman-simulation-nobias_tpm.rds')
saveRDS(t(sim_altman_noBias$cell_fractions), '/nfs/data/omnideconv_benchmarking_clean/data/PBMC/altman-simulation-nobias/altman-simulation-nobias_facs.rds')


#### Human brain ####

##### Myers #####
myers_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/matrix_counts.rds')
myers_cpm <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/matrix_norm_counts.rds')
cell_type <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/celltype_annotations.rds')

# SimBu dataset for simulations 
myers_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(myers_counts), cell_type = cell_type), 
  count_matrix = myers_counts,
  tpm_matrix = myers_cpm,
  name = 'Myers', 
  filter_genes = T)

saveRDS(myers_ds, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/simbu_ds.rds')

simulated_bulk <- SimBu::simulate_bulk(
  data =  myers_ds,
  scenario = 'random',
  scaling_factor = 'expressed_genes',
  nsamples = 50,
  ncells = 10000,
  total_read_counts = 100000000,  
  BPPARAM = BiocParallel::MulticoreParam(workers = 10),
  run_parallel = TRUE,  
  seed = 1234
)


##### Tran ######

tran_counts <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tran/matrix_counts.rds')
tran_cpm <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tran/matrix_norm_counts.rds')
cell_type <- readRDS('/nfs/data/omnideconv_benchmarking_clean/data/singleCell/tran/celltype_annotations.rds')

# SimBu dataset for simulations 
myers_ds <- SimBu::dataset(
  annotation = data.frame(ID = colnames(myers_counts), cell_type = cell_type), 
  count_matrix = myers_counts,
  tpm_matrix = myers_cpm,
  name = 'Myers', 
  filter_genes = T)

saveRDS(myers_ds, '/nfs/data/omnideconv_benchmarking_clean/data/singleCell/myers/simbu_ds.rds')




