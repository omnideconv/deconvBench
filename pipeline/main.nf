#!/usr/bin/env nextflow

include { PREPROCESS_SINGLE_CELL                            } from './modules/local/preprocess_single_cell'
include { SIMULATE_BULK                                     } from './modules/local/simulate_bulk'
include { CREATE_SIGNATURE                                  } from './modules/local/create_signature'
include { CREATE_SIGNATURE as CREATE_SIGNATURE_PREPROCESSED } from './modules/local/create_signature'
include { DECONVOLUTE                                       } from './modules/local/deconvolute'
include { COMPUTE_METRICS                                   } from './modules/local/compute_metrics'
include { ANALYSIS_BULK_MISSING_CELL_TYPES                   } from './modules/local/analysis_bulk_missing_cell_types'
include { SIMULATE_BULK_SPILLOVER                           } from './modules/local/simulate_bulk_spillover'
include { ANALYSIS_SPILLOVER                                } from './modules/local/analysis_spillover'
include { SIMULATE_BULK_UNKNOWN_CELL_TYPE                   } from './modules/local/simulate_bulk_unknown_cell_type'
include { ANALYSIS_BULK_UNKNOWN_CELL_TYPE                   } from './modules/local/analysis_bulk_unknown_cell_type'
include { SIMULATE_BULK_RESOLUTION_ANALYSIS                 } from './modules/local/simulate_bulk_resolution_analysis'
include { ANALYSIS_BULK_RESOLUTION_ANALYSIS                 } from './modules/local/analysis_bulk_resolution_analysis'
include { SIMULATE_PSEUDOBULK_MIRRORDB                      } from './modules/local/simulate_pseudobulk_mirrordb'
include { ANALYSIS_PSEUDOBULK_MIRRORDB                      } from './modules/local/analysis_pseudobulk_mirrordb'
include { ANALYSIS_BULK_MIRRORDB                            } from './modules/local/analysis_bulk_mirrordb'
include { SPLIT_REFERENCE                                   } from './modules/local/split_reference'
include { requireParams; datasetChannel; simulationInput; simulationDatasets; datasetsToSplit; deconvolutionInput; signatureInput; deconvolutionWithSignature; metricsInput } from './functions'

// The entry workflows call the processes directly instead of through sub-workflows, so that the
// process names in the trace files stay as they are (e.g. 'subsampling:DECONVOLUTE');
// visualisation/supplement/fig_s4.Rmd selects tasks by these names.

// Datasets for the simulation and for the methods. With split_reference, a dataset that has both
// roles is split by donor: one part goes to the simulation, the other to the methods. Without it,
// ch_to_split is empty and both channels pass through unchanged.
// All channels emit [name, directory].
workflow SPLIT_DATASETS {

  take:
  ch_simulation
  ch_reference
  ch_to_split

  main:
  SPLIT_REFERENCE(ch_to_split)

  emit:
  simulation = ch_simulation.join(SPLIT_REFERENCE.out.simulation, remainder: true).map { name, dir, part -> [name, part ?: dir] }
  reference  = ch_reference.join(SPLIT_REFERENCE.out.deconvolution, remainder: true).map { name, dir, part -> [name, part ?: dir] }
}

workflow {

  requireParams(['data_dir_sc', 'data_dir_bulk', 'results_dir_general', 'single_cell_list', 'bulk_list', 'method_list', 'ncores', 'species_sc'])

  ch_input = deconvolutionInput(datasetChannel(params.single_cell_list, params.data_dir_sc),
                                datasetChannel(params.bulk_list, params.data_dir_bulk),
                                false, true)

  CREATE_SIGNATURE(signatureInput(ch_input))

  DECONVOLUTE(deconvolutionWithSignature(ch_input, CREATE_SIGNATURE.out.signature))

  COMPUTE_METRICS(metricsInput(ch_input, DECONVOLUTE.out.results))
}

workflow deconv_only {

  requireParams(['data_dir_sc', 'data_dir_bulk', 'results_dir_general', 'single_cell_list', 'bulk_list', 'method_list', 'ncores', 'species_sc'])

  ch_input = deconvolutionInput(datasetChannel(params.single_cell_list, params.data_dir_sc),
                                datasetChannel(params.bulk_list, params.data_dir_bulk),
                                false, true)

  CREATE_SIGNATURE(signatureInput(ch_input))

  DECONVOLUTE(deconvolutionWithSignature(ch_input, CREATE_SIGNATURE.out.signature))
}

workflow simulation {

  requireParams(['data_dir_sc', 'data_dir_bulk', 'results_dir_general', 'single_cell_list', 'simulation_n_cells', 'simulation_n_samples', 'simulation_scenario', 'simulation_bias_type', 'method_list', 'ncores', 'species_sc'])

  simulation_names = simulationDatasets()
  SPLIT_DATASETS(datasetChannel(simulation_names, params.data_dir_sc),
                 datasetChannel(params.single_cell_list, params.data_dir_sc),
                 datasetChannel(datasetsToSplit(simulation_names, params.single_cell_list), params.data_dir_sc))

  SIMULATE_BULK(simulationInput(SPLIT_DATASETS.out.simulation, params.simulation_scenario, params.simulation_bias_type))

  ch_input = deconvolutionInput(SPLIT_DATASETS.out.reference,
                                SIMULATE_BULK.out.bulk,
                                false, true)

  CREATE_SIGNATURE(signatureInput(ch_input))

  DECONVOLUTE(deconvolutionWithSignature(ch_input, CREATE_SIGNATURE.out.signature))

  COMPUTE_METRICS(metricsInput(ch_input, DECONVOLUTE.out.results))
}

workflow subsampling {

  requireParams(['data_dir_sc', 'data_dir_bulk', 'results_dir_general', 'preProcess_dir', 'single_cell_list', 'ct_fractions', 'replicates', 'simulation_n_cells', 'simulation_n_samples', 'method_list', 'ncores', 'species_sc'])

  simulation_names = simulationDatasets()
  SPLIT_DATASETS(datasetChannel(simulation_names, params.data_dir_sc),
                 datasetChannel(params.single_cell_list, params.data_dir_sc),
                 datasetChannel(datasetsToSplit(simulation_names, params.single_cell_list), params.data_dir_sc))

  PREPROCESS_SINGLE_CELL(SPLIT_DATASETS.out.reference
                           .combine(channel.fromList(params.ct_fractions))
                           .combine(channel.of(1..params.replicates)))

  SIMULATE_BULK(simulationInput(SPLIT_DATASETS.out.simulation, ['random'], ['expressed_genes']))

  // the subsampled reference is named after its directory, e.g. sc1_perc100_rep1
  // Signatures are built for every bulk dataset and every method here, also where one signature
  // could be shared: fig_s4.Rmd shows the resources of the signature step per bulk dataset.
  ch_input = deconvolutionInput(PREPROCESS_SINGLE_CELL.out.subsample.map { subsample -> [subsample.name, subsample] },
                                SIMULATE_BULK.out.bulk,
                                true, false)

  CREATE_SIGNATURE_PREPROCESSED(signatureInput(ch_input))

  DECONVOLUTE(deconvolutionWithSignature(ch_input, CREATE_SIGNATURE_PREPROCESSED.out.signature))
}

workflow impact_missing_cell_types {

  requireParams(['data_dir_sc', 'data_dir_bulk', 'results_dir_missing_cell_types', 'single_cell_list', 'bulk_list', 'cell_types_to_exclude', 'method_list', 'ncores'])

  ch_bulk = datasetChannel(params.bulk_list, params.data_dir_bulk)
  ch_reference = datasetChannel(params.single_cell_list, params.data_dir_sc)

  // optional: also deconvolve the pseudo-bulks that SIMULATE_BULK creates from simulation_sc_dataset
  simulate = params.containsKey('simulation_missing_cell_types') && params.simulation_missing_cell_types
  simulation_names = []
  if (simulate) {
    requireParams(['preProcess_dir', 'simulation_n_cells', 'simulation_n_samples', 'simulation_scenario', 'simulation_bias_type'])
    simulation_names = simulationDatasets()
  }

  // With split_reference the references are split by donor in any case and the methods get the
  // method part, as in the workflows that simulate; without pseudo-bulks the other part is not used.
  split_names = simulate ? datasetsToSplit(simulation_names, params.single_cell_list) : datasetsToSplit(params.single_cell_list, params.single_cell_list)
  SPLIT_DATASETS(datasetChannel(simulation_names, params.data_dir_sc),
                 ch_reference,
                 datasetChannel(split_names, params.data_dir_sc))
  ch_reference = SPLIT_DATASETS.out.reference

  if (simulate) {
    SIMULATE_BULK(simulationInput(SPLIT_DATASETS.out.simulation, params.simulation_scenario, params.simulation_bias_type))
    ch_bulk = ch_bulk.mix(SIMULATE_BULK.out.bulk)
  }

  ANALYSIS_BULK_MISSING_CELL_TYPES(ch_reference
                                     .combine(ch_bulk)
                                     .combine(channel.fromList(params.method_list))
                                     .combine(channel.fromList(params.cell_types_to_exclude)))
}

workflow simulation_spillover {

  requireParams(['data_dir_sc', 'results_dir_spillover', 'preProcess_dir', 'single_cell_list', 'simulation_n_cells', 'simulation_n_samples', 'spillover_celltypes', 'method_list', 'ncores'])

  simulation_names = simulationDatasets()
  SPLIT_DATASETS(datasetChannel(simulation_names, params.data_dir_sc),
                 datasetChannel(params.single_cell_list, params.data_dir_sc),
                 datasetChannel(datasetsToSplit(simulation_names, params.single_cell_list), params.data_dir_sc))

  SIMULATE_BULK_SPILLOVER(SPLIT_DATASETS.out.simulation
                            .combine(channel.fromList(params.simulation_n_cells))
                            .combine(channel.fromList(params.simulation_n_samples)))

  ANALYSIS_SPILLOVER(SIMULATE_BULK_SPILLOVER.out.bulk
                       .combine(channel.fromList(params.method_list))
                       .combine(SPLIT_DATASETS.out.reference))
}

workflow simulation_unknown_content {

  requireParams(['data_dir_sc', 'results_dir_unknown_content', 'preProcess_dir', 'simulation_sc_dataset', 'simulation_n_cells', 'simulation_n_samples', 'known_cell_types', 'unknown_cell_type', 'fractions_unknown_cells', 'replicates_unknown_content', 'method_list', 'ncores'])

  // the dataset is used for the simulation and as reference
  ch_sc = datasetChannel(params.simulation_sc_dataset, params.data_dir_sc)
  SPLIT_DATASETS(ch_sc, ch_sc, datasetChannel(datasetsToSplit(params.simulation_sc_dataset, params.simulation_sc_dataset), params.data_dir_sc))

  SIMULATE_BULK_UNKNOWN_CELL_TYPE(SPLIT_DATASETS.out.simulation
                                    .combine(channel.fromList(params.simulation_n_cells))
                                    .combine(channel.fromList(params.simulation_n_samples))
                                    .combine(channel.fromList(params.unknown_cell_type))
                                    .combine(channel.fromList(params.replicates_unknown_content)))

  // the pseudo-bulks are deconvolved with the dataset they were simulated from
  ANALYSIS_BULK_UNKNOWN_CELL_TYPE(SIMULATE_BULK_UNKNOWN_CELL_TYPE.out.bulk
                                    .combine(SPLIT_DATASETS.out.reference, by: 0)
                                    .combine(channel.fromList(params.method_list)))
}

workflow simulation_resolution_analysis {

  requireParams(['data_dir_sc', 'results_dir_resolution', 'preProcess_dir', 'simulation_sc_dataset', 'simulation_n_cells', 'simulation_n_samples', 'cell_types_finer_res', 'replicates_simulation', 'method_list', 'ncores'])

  // the dataset is used for the simulation and as reference
  ch_sc = datasetChannel(params.simulation_sc_dataset, params.data_dir_sc)
  SPLIT_DATASETS(ch_sc, ch_sc, datasetChannel(datasetsToSplit(params.simulation_sc_dataset, params.simulation_sc_dataset), params.data_dir_sc))

  SIMULATE_BULK_RESOLUTION_ANALYSIS(SPLIT_DATASETS.out.simulation
                                      .combine(channel.fromList(params.simulation_n_cells))
                                      .combine(channel.fromList(params.simulation_n_samples))
                                      .combine(channel.fromList(params.cell_types_finer_res))
                                      .combine(channel.fromList(params.replicates_simulation)))

  // the pseudo-bulks are deconvolved with the dataset they were simulated from
  ANALYSIS_BULK_RESOLUTION_ANALYSIS(SIMULATE_BULK_RESOLUTION_ANALYSIS.out.bulk
                                      .combine(SPLIT_DATASETS.out.reference, by: 0)
                                      .combine(channel.fromList(params.method_list)))
}

workflow simulation_impact_technology {

  requireParams(['data_dir_sc', 'data_dir_bulk', 'results_dir_impact_technology', 'preProcess_dir', 'simulation_sc_dataset', 'simulation_n_cells', 'simulation_n_samples', 'sc_datasets_impact_technology', 'bulk_impact_technology', 'replicates_simulation', 'method_list', 'ncores'])

  // sc_datasets_impact_technology is a list with one entry of comma-separated dataset names
  all_names = params.sc_datasets_impact_technology.collectMany { names -> names.tokenize(',') }

  // pseudo-bulks are simulated from all datasets; simulation_sc_dataset are the references
  SPLIT_DATASETS(datasetChannel(all_names, params.data_dir_sc),
                 datasetChannel(params.simulation_sc_dataset, params.data_dir_sc),
                 datasetChannel(datasetsToSplit(all_names, params.simulation_sc_dataset), params.data_dir_sc))

  ch_sc = SPLIT_DATASETS.out.reference
  ch_all_sc_dirs = SPLIT_DATASETS.out.simulation
                     .map { _name, sc_dir -> sc_dir }
                     .collect()

  SIMULATE_PSEUDOBULK_MIRRORDB(channel.fromList(params.simulation_sc_dataset)
                                 .combine(channel.fromList(params.simulation_n_cells))
                                 .combine(channel.fromList(params.simulation_n_samples))
                                 .combine(channel.fromList(params.replicates_simulation)),
                               ch_all_sc_dirs)

  // every reference is used on the pseudo-bulks of every dataset
  ANALYSIS_PSEUDOBULK_MIRRORDB(ch_sc
                                 .combine(SIMULATE_PSEUDOBULK_MIRRORDB.out.bulk)
                                 .combine(channel.fromList(params.method_list)),
                               SIMULATE_PSEUDOBULK_MIRRORDB.out.cell_types.first())

  ANALYSIS_BULK_MIRRORDB(ch_sc
                           .combine(datasetChannel(params.bulk_impact_technology, params.data_dir_bulk))
                           .combine(channel.fromList(params.method_list)))
}
