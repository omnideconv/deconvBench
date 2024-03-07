# deconvBench: a nextflow pipeline to benchmark second-generation deconvolution methods using omnideconv

For our effort to benchmark eight second-generation deconvolution methods, all included in [omnideconv](https://github.com/omnideconv/omnideconv), we implemented a nextflow pipeline that includes several workflows. These offer various simulation scenarios and the option to deconvolve real bulk RNA-seq datasets, all with different scRNA-seq datasets, deconvolution methods and other parameters. 
A detailed description of the parameters and workflows follows below.

## Dependencies

All method-related dependencies are installed in a inlcuded Docker image, that you can either build yourself using the Dockerfile in the directory above (`docker/`) or pull it directly from [Dockerhub](https://hub.docker.com/repository/docker/alexd13/omnideconv_benchmark/general).

In addition, to run the pipeline you need an installation of [nextflow](https://www.nextflow.io/). 

## Quickstart

If you want to deconvolve two bulk RNA-seq datasets (`bulk1` and `bulk2`) using one scRNA-seq dataset (`sc1`) and compute quality metrics for four selected deconvolution methods (`AutoGeneS`,`DWLS`,`MuSiC` and `Scaden`), you can start the pipeline like this, inside the current directory:

```
nextflow -C nextflow.config run main.nf -profile docker
```

The `nextflow.config` file can then look something like this:

```
params {
    data_dir_bulk = "/user/benchmarking/datasets/bulks"
    data_dir_sc = "/user/benchmarking/datasets/single_cell"
    results_dir_general = "/user/benchmarking/results"   

    single_cell_list = ["sc1"]
    bulk_list = ["bulk1","bulk2"]
    method_list = ["autogenes", "dwls", "music", "scaden"]

    ncores = '4'
}

profiles {
    docker {
        docker.enabled = true
        process.container = 'omnideconv_benchmark:latest'
        // add more parameters regarding e.g. memory usage here
    }
}

```
Running the pipeline with the above configuration will expect the following folder structure:


```
user
└── benchmarking
    ├── datasets
    │   ├── bulks
    │   │   ├── bulk1
    |   |   |   ├── bulk1_counts.rds
    |   |   |   ├── bulk1_tpm.rds
    |   |   |   └── bulk1_facs.rds
    │   │   └── bulk2
    |   |       └── ...
    │   └── single_cell
    │       └── sc1
    │           ├── batch.rds
    │           ├── celltype_annotation.rds
    │           ├── matrix_counts.rds
    │           └── matrix_norm_counts.rds
    ├── results
    │
    └── DECONVBENCH
        ├── docker
        ├── pipeline
        │   ├── nextflow.config
        │   ├── main.nf
        │   ├── optimal_normalization.csv
        │   ├── bin
        │   │   ├── runDeconvolutionNF.R
        │   │   ├── ...
        └── visualisation
```


## Workflows

![Alt text](workflows.png)

### main
The **main** workflow is run by using nextflow without the `-entry` parameter.


This workflow does signature creation, deconvolution and result metric calculations on the selected deconvolution methods (`method_list`) and datasets (`single_cell_list`, `bulk_list`).

In the **main** workflow following processes are run:
1. CREATE_SIGNATURE - Creates signature matrix from scRNA-seq datasets set in `single_cell_list`.
2. DECONVOLUTION - Deconvolves the bulk RNA-seq datasets set in `bulk_list` using the  created signature matrices.
3. COMPUTE_METRICS - Calculates the metrics for the deconvolution results.

The results of the **main** workflow are stored in the `results_dir_general` directory.

### subsampling
The **subsampling** workflow is run by using nextflow with the `-entry subsampling` parameter.


This workflow builds upon the **main** workflow and adds a cell-type specific subsampling step. The size of the subsets and the number of replicates can be controlled by the parameters `ct_fractions` and `replicates`.

In the **subsampling** workflow following processes are run:
1. PREPROCESS_SINGLE_CELL - For each scRNA-seq dataset in `single_cell_list` the dataset is subsampled and preprocessed. If `ct_fractions` contains integer, the number of cells to be subsampled per celltype is equal to the integer, which corresponds to an *even* sampling. For floats, we do a *mirrordb* style sampling, where  a corresponding fraction of the total number of cells per celltype is sampled. The number of replicates is controlled by the `replicates` parameter. Cells are drawn randomly without replacement.
2. CREATE_SIGNATURE - Creates signature matrix from the subsampled scRNA-seq datasets.
3. DECONVOLUTION - Deconvolves the bulk RNA-seq datasets set in `bulk_list` using the  created signature matrices.
4. COMPUTE_METRICS - Calculates the metrics for the deconvolution results.

The results of the **subsampling** workflow are stored in the `results_dir_general` directory.

### impact_missing_cell_types
The **impact_missing_cell_types** workflow is run by using nextflow with the `-entry impact_missing_cell_types` parameter.


This workflow excludes cell types from the scRNA-seq datasets (`single_cell_list`) and then deconvolves the bulk RNA-seq datasets(`bulk_list`) using the created signature matrices. The cell types to be excluded are controlled by the `cell_types_to_exclude` parameter. The workflow is repeated for each cell type in the list, where each time one cell type is excluded.

In the **impact_missing_cell_types** workflow following processes is run:
1. ANALYSIS_BULK_MISSING_CELL_TYPES - Excludes `cell_types_to_exclude` cell types from the `single_cell_list` scRNA-seq datasets and deconvolves the bulk RNA-seq datasets using the created signature matrices. Finally, the metrics for the deconvolution results are calculated.

The results of the **impact_missing_cell_types** workflow are stored in the `results_dir_missing_cell_types` directory.

### simulation_unknown_content

The **simulation_unknown_content** workflow is run by using nextflow with the `-entry simulation_unknown_content` parameter.


This workflow simulates the presence of an unknown cell type in  the bulk RNA-seq datasets (`bulk_list`). The unknown cell type is controlled by the `unknown_cell_type` parameter.

In the **simulation_unknown_content** workflow following processes are run:
1. SIMULATE_BULK_UNKNOWN_CELL_TYPE - Creates simulated bulk RNA-seq datasets with the presence of `unknown_cell_type` . The datasets are created using the `bulk_list` datasets. The fraction of the unknown cell type is controlled by the `fractions_unknown_cells` parameter. The number of replicates is controlled by the `replicates_unknown_content` parameter. 

2. ANALYSIS_BULK_UNKNOWN_CELL_TYPE - Creates signature matrices using `known_cell_types` subsetss of the `single_cell_list` datasets and deconvolves the simulated bulk RNA-seq datasets using the created signature matrices. Finally, the metrics for the deconvolution results are calculated.

The results of the **simulation_unknown_content** workflow are stored in the `results_dir_unknown_content` directory.


### simulation_impact_technology

The **simulation_impact_technology** workflow is run by using nextflow with the `-entry simulation_impact_technology` parameter.

In this workflow we evaluate the impact of differences in tissue origin and sequencing technology on deconvolution results. 

In the **simulation_impact_technology** workflow following processes are run:
1. SIMULATE_PSEUDOBULK_MIRRORDB - 

2. ANALYSIS_PSEUDOBULK_MIRRORDB - 

The results of the **simulation_impact_technology** workflow are stored in the `results_dir_impact_technology` directory.

### simulation_resolution_analysis

The **simulation_resolution_analysis** workflow is run by using nextflow with the `-entry simulation_resolution_analysis` parameter.

This workflow assesses the impact of annoating the usage of single cell atlases with different levels of cell type annotations. 

In the **simulation_resolution_analysis** workflow following processes are run:

1. SIMULATE_BULK_RESOLUTION_ANALYSIS -

2. ANALYSIS_BULK_RESOLUTION_ANALYSIS -


The results of the **simulation_resolution_analysis** workflow are stored in the `results_dir_resolution` directory.


### simulation_spillover


The **simulation_spillover** workflow is run by using nextflow with the `-entry simulation_spillover` parameter.

This workflow assesses the impact  the erroneous prediction of a cell type caused by the occurance of cell types with similar transcriptional profiles.

In the **simulation_spillover** workflow following processes are run:

1. SIMULATE_BULK_SPILLOVER -

2. ANALYSIS_SPILLOVER -

The results of the **simulation_spillover** workflow are stored in the `results_dir_spillover` directory.


## Parameters

Parameters can be changed in the `nextflow.config` file.

### Input Directories

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `data_dir_bulk` | path | absolute path to directory that contains RNA-seq datasets | `/path/to/datasets/` |
| `data_dir_sc` | path | absolute path to directory that contains scRNA-seq datasets                          | `/path/to/other/datasets/`          |

A dedicted description of the input file formats and folder structure can be found [below](#input-formats).

### Output Directories

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `results_dir_general` | path | absolute path to directory in which final results of the `main` and `subsampling` workflow are stored | `/path/to/results/` |
| `preProcess_dir` | path | absolute path to directory in which temporary files are stored                          | `/path/to/temporary/directory/`          |

### General Parameters

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `ncores` | int | Number of cores that are available for methods | `4` |
| `species_sc` | string | Type of species in scRNA-seq dataset (`mm` or `hs`). Currently only supported by BayesPrism. | `hs` |

### Benchmarking Parameters [general]

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `single_cell_list` | [string] | list of scRNA-seq dataset names that will be used for deconvolution. Details for naming conventions are found [below](#input-formats). | `["hao-complete", "hao-sampled-3"]` |
| `bulk_list` | [string] |   list of RNA-seq dataset names that will be used for deconvolution. Details for naming conventions are found [below](#input-formats).                        |    `["finotello", "finotello-simulation"]`      |
| `method_list` | [string] |  list of method names that are used for deconvolution. Possible values are: `["autogenes","bayesprism","bisque","cibersortx"` `,"dwls","music","scaden","scdc"]`                         |    `["music","bayesprism"]`       |

### Workflow-specific Parameters:

The following sets of parameters are only used in the specified workflow.

####  workflow `subsampling` 

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `ct_fractions` | [numeric] | This parameter influences how large the subsets are that are drawn from the scRNA-seq dataset. There are two ways of using it: (1) Absolute number of cells, indicated by integers `>= 1`. If there are less cells available for one cell type than the number suggests, all cells of this type are selected. Otherwise, cells are drawn randomly without replacement. (2) Relative fraction of cells, indicated by numeric values in the range of `[0.001-0.999]`. This selects a subset of cells for each cell type, relative to the total number of cells in this type. | `[0.01, 0.05, 100, 500]` |
| `replicates` | int | The number of technical replicates that are drawn from the scRNA-seq dataset. Repeats the subsampling with the indicated `ct_fractions` the indicated amount of times, in order to account for technical biases during random cell sampling. | `5` |


#### workflow `simulation_impact_technology` 

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `results_dir_impact_technology` | path | absolute path to directory in which final results of the `impact_technology` workflow are stored | `/path/to/impact_technology/results/` |
| `replicates_simulation` | [int] | number of technical replicates for pseudo-bulk simulation | `[5]` |
| `datasets_impact_technology` | [string] | names of scRNA-seq datasets in `data_dir_sc` that are used for this workflow | `["lambrechts","maynard]` |


#### workflow `simulation_spillover` 

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `results_dir_spillover` | path | absolute path to directory in which final results of the `simulation_spillover` workflow are stored | `/path/to/simulation_spillover/results/` |
| `spillover_samples_per_cell` | int | in the spillover analysis samples contain only one cell type. This parameter controls how many pseudo-bulk samples should be simulated that way. | `10` |
| `spillover_celltypes` | [string] | Control, for which cell types the spillover analysis will be done | `["B cells,Monocytes,NK cells,T cells CD8,T cells CD4 conv"]` |

#### workflow `simulation_unkown_content` 

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `results_dir_unknown_content` | path | absolute path to directory in which final results of the `simulation_unkown_content` workflow are stored | `/path/to/simulation_unkown_content/results/` |
| `known_cell_types` | [[string]] | subset of cell types that we will use to build the signature matrix | `[["B cells,Stromal cells,T cells CD4 conv,Macrophages"]]` |
| `unknown_cell_type` | [string] | unknown cell type to use | `["Tumor cells"]` |
| `fractions_unknown_cells` | [[int]] | Control, for which cell types the spillover analysis will be done | `[[0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9]]` |
| `replicates_unknown_content` | [int] | number of technical replicates for pseudo-bulk simulation | `[5]` |

#### workflow `simulation_resolution_analysis` 

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `results_dir_resolution` | path | absolute path to directory in which final results of the `simulation_resolution_analysis` workflow are stored | `/path/to/simulation_resolution_analysis/results/` |
| `cell_types_finer_res` | [[string]] | which cell types will we inspect at various resolutions? | `["Endothelial ACKR1,Endothelial RGS5,Endothelial CXCL12,CAFs MSC iCAF-like s1,CAFs MSC iCAF-like s2, ..."]` |

#### workflow `impact_missing_cell_types` 

| Name       | Type   | Description                             | Example                      |
| ---------- | ------ | ---------------------------------- | -------------------------------- |
| `results_dir_missing_cell_types` | path | absolute path to directory in which final results of the `impact_missing_cell_types` workflow are stored | `/path/to/impact_missing_cell_types/results/` |
| `cell_types_to_exclude` | [string] | each of the cell types in the list will be a candidate for exclusion from the scRNA-seq dataset once | `["B cells","mDC","Monocytes","NK cells","T cells CD4 conv"]` |

## Input Formats

## Data Normalizations

Each deconvolution method has a recommended set of normalization procedures, both for the scRNA-seq and bulk RNA-seq samples. Possible options currently are `counts` (un-normalized), `tpm` and `cpm`. The specific values for each method are stored in the `optimal_normalization.csv` file, and at time of publication look like this:

| method       | sc_norm   | bulk_norm |
|----|----|----|
|autogenes|cpm|tpm|
|bayesprism|counts|counts|
|bisque|counts|counts|
|cibersortx|cpm|tpm|
|dwls|counts|tpm|
|music|counts|tpm|
|scaden|counts|tpm|
|scdc|counts|tpm|
