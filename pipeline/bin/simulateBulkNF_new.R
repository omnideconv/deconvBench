#!/usr/bin/Rscript
print("start simulator")

"Usage: 
  simulateBulkNF.R <sc_path> <sc_meta> <remapping_sheet> <scenario> [<number> <tissue> <organism> <celltypes>]
Options:
<sc_path> path to sc dataset
<sc_meta> metainformation and paths to sc data --> [[maynard_tpm:maynard/matrix_tpm.rds, maynard_counts:maynard/matrix_counts.rds], [lambrechts_tpm:/lambrechts/matrix_norm_counts.rds, lambrechts_counts:lambrechts/matrix_counts.rds]]
<scenario> simulation scenario (spillover, uniform, etc)
<number> spike in percent
<celltypes> if spillover is chosen, give list of celltypes (correct names)" -> doc

args <- docopt::docopt(doc)
number <- as.integer(args$number)

sc_path <- args$sc_path
sc_meta <- strsplit(gsub("\\]", "", gsub("\\[", "", args$sc_meta)), split = ", ")[[1]]
sc_ds <- strsplit(sc_meta[1], split="_")[[1]][1]
meta_list <- list()
for(i in 1:length(sc_meta)){
  split <- strsplit(sc_meta[i], split=":")[[1]];
  names_ml <- names(meta_list)
  meta_list <- c(meta_list, split[2])
  names(meta_list) <- c(names_ml, strsplit(split[1], "_")[[1]][2])
}
sc_matrix_tpm <- readRDS(meta_list$tpm)
sc_matrix_counts <- readRDS(meta_list$counts)
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))
scenario <- args$scenario

remapping_sheet <- args$remapping_sheet

tissue <- args$tissue
organism <- args$organism

source("/vol/spool/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)
celltypes <- unlist(strsplit(args$celltypes, ","))
celltypes <- unique(sc_celltype_annotations)#[1:3] #remove this if you only want certain celltypes for spillover!
#celltypes <- c("B cell")

escapeCelltypes <- function(celltype){
  return(gsub(" ", "x.x", celltype))
}
reEscapeCelltypes <- function(celltype){
  return(gsub("x.x", " ", celltype))
}

library(SimBu)
library(data.table)
library(Matrix)
library(BiocParallel)
if(scenario=="uniform"){
  ###even###
  simulated_bulk <- simulate_bulk(dataset(
    annotation = data.frame(ID = colnames(sc_matrix_counts), cell_type = sc_celltype_annotations), 
    count_matrix = sc_matrix_counts,
    tpm_matrix = sc_matrix_tpm,
    name = sc_ds), scenario = "even")
} else if(scenario=="truefractions"){
  ###mirror_db###
  simulated_bulk <- simulate_bulk(dataset(
    data.frame(ID = colnames(sc_matrix_counts), cell_type = sc_celltype_annotations), 
    count_matrix = sc_matrix_counts,
    tpm_matrix = sc_matrix_tpm,
    name = sc_ds), scenario = "mirror_db")
} else if(scenario == "spillover"){
  ###spillover###
  for (celltype in celltypes) {
    print(celltype)
    simulated_bulk <- SimBu::simulate_bulk(
      dataset(
      annotation = data.frame(ID = colnames(sc_matrix_counts), 
                              cell_type = sc_celltype_annotations),
      count_matrix = sc_matrix_tpm,
      tpm_matrix = sc_matrix_tpm,
      name = sc_ds),
      scaling_factor = "NONE",
      scenario = "pure", 
      pure_cell_type=celltype)
    print("save simulation")
    saveRDS(simulated_bulk, 
            paste("simulatedBulk_", scenario, "_", sc_ds, "_", escapeCelltypes(celltype), ".rds", sep=""), 
            compress = FALSE)
  }
  quit()
} else if(scenario == "spikein"){
  ###spikein###
  celltype="T cell CD8+"
  #for (celltype in celltypes) {
    print(celltype)
    #first create the background
    simulated_bulk <- SimBu::simulate_bulk(
      dataset(annotation = data.frame(ID = colnames(sc_matrix_counts), 
                                      cell_type = sc_celltype_annotations),
        count_matrix = sc_matrix_tpm,
        tpm_matrix = sc_matrix_tpm,
        name = sc_ds),
      scaling_factor = "NONE",
      scenario = "random", 
      nsamples = 1000,
      blacklist = c(celltype), BPPARAM=bpparam("MulticoreParam"), run_parallel = TRUE)
    print("save simulation")
    saveRDS(simulated_bulk, 
            paste("simulatedBulk_", scenario, "_", sc_ds, "_backgroundFor_", escapeCelltypes(celltype), ".rds", sep=""), 
            compress = FALSE)
    #make spike in sets with increasing amount of that cell (param number)
    #for(number in c(1, 2, 3, 5, 10, 15, 20, 30, 50, 100)){ #, 20, 30, 50
      for(x in 1:5){
        simulated_bulk <- SimBu::simulate_bulk(
          dataset(
            annotation = data.frame(ID = colnames(sc_matrix_counts), 
                                    cell_type = sc_celltype_annotations),
            count_matrix = sc_matrix_tpm,
            tpm_matrix = sc_matrix_tpm,
            name = sc_ds),
          scaling_factor = "NONE",
          scenario = "weighted", 
          nsamples = 1000,
          weighted_cell_type = celltype, 
          weighted_amount = number/1000, 
          BPPARAM=bpparam("MulticoreParam"), run_parallel = TRUE)
        print("save simulation")
        saveRDS(simulated_bulk, 
                paste("simulatedBulk_", scenario, "_", sc_ds, "_spikeinFor_", escapeCelltypes(celltype), "_nrCells", number, "_x", x,".rds", sep=""), 
                compress = FALSE)
      }
    #}
  #}
  quit()
} else if(scenario=="sfaira"){
  setup_list <- SimBu::setup_sfaira(
    python_path = "/nfs/home/students/k.reinisch/.conda/envs/sfaira/bin/python3",
    env_name = "sfaira", 
    basedir = "/nfs/data/omnideconv_benchmarking/data/sfaira/")
  ds_human_blood <- SimBu::dataset_sfaira_multiple(sfaira_setup = setup_list,
                                                  organisms = "Homo sapiens", 
                                                  tissues = "blood", 
                                                  assays = "CITE-seq (cell surface protein profiling)",
                                                  name="human_blood")
  ds_single <- SimBu::dataset_sfaira(sfaira_setup = setup_list, 
                                     sfaira_id = "homosapiens_blood_2020_citeseq(cellsurfaceproteinprofiling)_haoyuhan_001_10.1016/j.cell.2021.04.048")
  
  simulated_bulk <- simulate_bulk(dataset(
    annotation = data.frame(ID = colnames(sc_matrix_counts), cell_type = sc_celltype_annotations),
    count_matrix = sc_matrix_counts,
    tpm_matrix = sc_matrix_tpm,
    name = sc_ds),
    scaling_factor = "NONE",
    scenario = "unique", 
    unique_cell_type="B cell")
}

saveRDS(simulated_bulk, 
        paste("simulatedBulk_", scenario, "_", sc_ds, ".rds", sep=""), 
        compress = FALSE)

