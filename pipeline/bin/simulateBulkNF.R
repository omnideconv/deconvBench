#!/usr/bin/Rscript
print("start simulator")

"Usage: 
  simulateBulkNF.R <sc_path> <sc_datasetname> <remapping_sheet> <scenario> [<tissue> <organism> <celltypes>]
Options:
<sc_path> path to sc dataset
<sc_datasetname> name of sc dataset
<remapping_sheet> excel mapping sheet with remapping
<scenario> simulation scenario (spillover, uniform, sfaira)
<tissue> if sfaira is chosen, you can choose the tissue
<organism> if sfaira is chosen, you can choose the organism
<celltypes> if spillover is chosen, give list of celltypes (correct names)" -> doc

args <- docopt::docopt(doc)

sc_path <- args$sc_path
sc_ds <- args$sc_datasetname
sc_matrix <- readRDS(file.path(sc_path, sc_ds, "matrix.rds"))
sc_celltype_annotations <- readRDS(file.path(sc_path, sc_ds, "celltype_annotations.rds"))
sc_batch <- readRDS(file.path(sc_path, sc_ds, "batch.rds"))
sc_marker <- readRDS(file.path(sc_path, sc_ds, "marker.rds"))
scenario <- args$scenario

remapping_sheet <- args$remapping_sheet

tissue <- args$tissue
organism <- args$organism
celltypes <- unlist(strsplit(args$celltypes, ","))
celltypes <- unique(sc_celltype_annotations) #remove this if you only want certain celltypes for spillover!

source("/nfs/proj/omnideconv_benchmarking/benchmark/pipeline/bin/remapCelltypesNF.R")
sc_celltype_annotations <- remapCelltypesWorkflow(remappingPath = remapping_sheet, 
                                                  celltype_annotations = sc_celltype_annotations, 
                                                  method_ds = sc_ds)

library(SimBu)
library(data.table)
if(scenario=="uniform"){
  ###uniform###
  simulated_bulk <- simulate_bulk(dataset(
    annotation = data.frame(ID = colnames(sc_matrix), cell_type = sc_celltype_annotations), 
    count_matrix = sc_matrix,
    name = sc_ds))
} else if(scenario=="true_fractions"){
  ###mirror_db###
  simulated_bulk <- simulate_bulk(dataset(
    data.frame(ID = colnames(sc_matrix), cell_type = sc_celltype_annotations), 
    count_matrix = sc_matrix,
    name = sc_ds), scenario = "mirror_db")
} else if(scenario == "spillover"){
  ###spillover###
  for (celltype in celltypes) {
    print(celltype)
    simulated_bulk <- SimBu::simulate_bulk(
      dataset(
      annotation = data.frame(ID = colnames(sc_matrix), 
                              cell_type = sc_celltype_annotations),
      count_matrix = sc_matrix,
      name = sc_ds),
      scaling_factor = "NONE",
      scenario = "unique", 
      unique_cell_type=celltype) #TODO make this variable
    saveRDS(simulated_bulk, 
            paste("simulatedBulk_", scenario, "_", sc_ds, "_", celltype, ".rds", sep=""), 
            compress = FALSE)
  }
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
    annotation = data.frame(ID = colnames(sc_matrix), cell_type = sc_celltype_annotations),
    count_matrix = sc_matrix,
    name = sc_ds),
    scaling_factor = "NONE",
    scenario = "unique", 
    unique_cell_type="B cell")
}

saveRDS(simulated_bulk, 
        paste("simulatedBulk_", scenario, "_", sc_ds, ".rds", sep=""), 
        compress = FALSE)

