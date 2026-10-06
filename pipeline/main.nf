#!/usr/bin/env nextflow

include { PREPROCESS_SINGLE_CELL                            } from './modules/local/preprocess_single_cell'
include { SIMULATE_BULK                                     } from './modules/local/simulate_bulk'
include { CREATE_SIGNATURE                                  } from './modules/local/create_signature'
include { CREATE_SIGNATURE as CREATE_SIGNATURE_PREPROCESSED } from './modules/local/create_signature'
include { DECONVOLUTE                                       } from './modules/local/deconvolute'
include { COMPUTE_METRICS                                   } from './modules/local/compute_metrics'
include { datasetChannel; simulationInput; deconvolutionInput; metricsInput } from './functions'

// The entry workflows call the processes directly instead of through sub-workflows, so that the
// process names in the trace files stay as they are (e.g. 'subsampling:DECONVOLUTE');
// visualisation/supplement/fig_s4.Rmd selects tasks by these names.

/*
 * Workflows that pass files between processes
 */

workflow {

  ch_input = deconvolutionInput(datasetChannel(params.single_cell_list, params.data_dir_sc),
                                datasetChannel(params.bulk_list, params.data_dir_bulk),
                                false)

  CREATE_SIGNATURE(ch_input)

  DECONVOLUTE(ch_input.join(CREATE_SIGNATURE.out.signature))

  COMPUTE_METRICS(metricsInput(ch_input, CREATE_SIGNATURE.out.signature, DECONVOLUTE.out.deconvolution))
}

workflow deconv_only {

  ch_input = deconvolutionInput(datasetChannel(params.single_cell_list, params.data_dir_sc),
                                datasetChannel(params.bulk_list, params.data_dir_bulk),
                                false)

  CREATE_SIGNATURE(ch_input)

  DECONVOLUTE(ch_input.join(CREATE_SIGNATURE.out.signature))
}

workflow simulation {

  SIMULATE_BULK(simulationInput(params.simulation_scenario, params.simulation_bias_type))

  ch_input = deconvolutionInput(datasetChannel(params.single_cell_list, params.data_dir_sc),
                                SIMULATE_BULK.out.bulk,
                                false)

  CREATE_SIGNATURE(ch_input)

  DECONVOLUTE(ch_input.join(CREATE_SIGNATURE.out.signature))

  COMPUTE_METRICS(metricsInput(ch_input, CREATE_SIGNATURE.out.signature, DECONVOLUTE.out.deconvolution))
}

workflow subsampling {

  PREPROCESS_SINGLE_CELL(datasetChannel(params.single_cell_list, params.data_dir_sc)
                           .combine(channel.fromList(params.ct_fractions))
                           .combine(channel.of(1..params.replicates)))

  SIMULATE_BULK(simulationInput(['random'], ['expressed_genes']))

  // the subsampled reference is named after its directory, e.g. sc1_perc100_rep1
  ch_input = deconvolutionInput(PREPROCESS_SINGLE_CELL.out.subsample.map { subsample -> [subsample.name, subsample] },
                                SIMULATE_BULK.out.bulk,
                                true)

  CREATE_SIGNATURE_PREPROCESSED(ch_input)

  DECONVOLUTE(ch_input.join(CREATE_SIGNATURE_PREPROCESSED.out.signature))
}

/*
 * Workflows that still pass dataset names; their scripts write directly into the result directories
 */

process ANALYSIS_BULK_MISSING_CELL_TYPES {

      label 'process_default'

      input:
      each sc_dataset
      each bulk_ds
      each method 
      each cell_types_missing
           
      output:
      val("${sc_dataset}-missing_cell_types")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_missing_cell_type.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{bulk_ds}' '!{params.data_dir_bulk}' '!{method}' '!{cell_types_missing}' '!{params.results_dir_missing_cell_types}' '!{params.ncores}' '!{baseDir}'
      ''' 
}

process SIMULATE_BULK_SPILLOVER {

      label 'process_default'
      
      publishDir "${params.preProcess_dir}/pseudo_bulk_spillover", mode: 'copy'

      input:
      each sc_dataset
      each simulation_n_cells
      each simulation_n_samples
      val cell_types

      output:
      tuple val("${sc_dataset}_spillover_sim"),
            val("${params.preProcess_dir}/pseudo_bulk_spillover/${sc_dataset}_spillover_sim")

      beforeScript 'chmod o+rw .'      
      shell:
      '''
      simulateBulkNF_spillover_analysis.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{cell_types}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_SPILLOVER {

      label 'process_high'

      input:
      tuple val(sim_bulk_name),
            val(sim_bulk_path)
      val cell_types
      each method 
      each sc_reference
     
      output:
      val("${sc_dataset}-spillover")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_spillover.R '!{sc_reference}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{method}' '!{cell_types}' '!{params.results_dir_spillover}' '!{params.ncores}' '!{baseDir}'
      ''' 
}

process SIMULATE_BULK_UNKNOWN_CELL_TYPE {
      
      publishDir "${params.preProcess_dir}/pseudo_bulk_unknown_content", mode: 'copy'

      input:
      each sc_dataset
      each simulation_n_cells
      each simulation_n_samples  
      val(fraction_unknown_cell)
      val(cell_types)
      each unknown_cell_type
      each replicates

      output:
      tuple val(sc_dataset),
            val("${sc_dataset}_unknown_content_sim"),
            val("${params.preProcess_dir}/pseudo_bulk_unknown_content/${sc_dataset}_unknown_content_sim")
            val(replicates)

      beforeScript 'chmod o+rw .'      
      shell:
      '''
      simulateBulkNF_unknown_content_analysis.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{fraction_unknown_cell}' '!{cell_types}' '!{unknown_cell_type}' '!{replicates}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_BULK_UNKNOWN_CELL_TYPE {

      input:
      tuple val(sc_dataset), 
            val(sim_bulk_name), 
            val(sim_bulk_path)
            val(replicates)
      val fraction_unknown_cell
      val cell_types
      each unknown_cell_type
      each method
     
      output:
      val("${sc_dataset}-spillover")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_unknown_content.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{replicates}' '!{method}' '!{fraction_unknown_cell}' '!{cell_types}' '!{unknown_cell_type}' '!{params.results_dir_unknown_content}' '!{params.ncores}' '!{baseDir}'
      ''' 
}

process SIMULATE_BULK_RESOLUTION_ANALYSIS {

      label 'process_default'
  
      publishDir "${params.preProcess_dir}/pseudo_bulk_resolution", mode: 'copy'

      input:
      each sc_dataset
      each simulation_n_cells
      each simulation_n_samples                   
      each cell_types_fine
      each replicates
      output:
      tuple val(sc_dataset),
            val("${sc_dataset}_resolution_analysis_sim"),
            val("${params.preProcess_dir}/pseudo_bulk_resolution/${sc_dataset}_resolution_analysis_sim"),
            val(replicates)
      
      beforeScript 'chmod o+rw .'     

      shell:
      '''
      simulateBulkNF_impact_cell_resolution.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{cell_types_fine}' '!{replicates}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_BULK_RESOLUTION_ANALYSIS {
      
      label 'process_high'

      publishDir "${params.preProcess_dir}/pseudo_bulk_resolution", mode: 'copy'

      input:
      tuple val(sc_dataset), 
            val(sim_bulk_name), 
            val(sim_bulk_path),
            val(replicates)
      val cell_types_fine
      each method
     
      output:
      val("${params.simulation_sc_dataset}-resolution_analysis")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_impact_cell_resolution.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{method}' '!{cell_types_fine}' '!{replicates}' '!{params.results_dir_resolution}' '!{params.ncores}' '!{baseDir}'
      ''' 
}

process SIMULATE_PSEUDOBULK_MIRRORDB {

      label 'process_high'
      
      publishDir "${params.preProcess_dir}/pseudo_bulk_impact_technology", mode: 'copy'

      input:
      each sc_dataset
      each simulation_n_cells
      each simulation_n_samples  
      val(sc_dataset_list)
      each replicates

      output:
      tuple val(sc_dataset),
            val("${params.preProcess_dir}/pseudo_bulk_impact_technology/${sc_dataset}"),
            val(replicates)

      beforeScript 'chmod o+rw .'      
      shell:
      '''
      simulateBulkNF_mirrordb.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{sc_dataset_list}' '!{replicates}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_PSEUDOBULK_MIRRORDB {

      label 'process_default'

      input:
      each sc_dataset
      tuple val(sim_bulk_name), 
            val(sim_bulk_path),
            val(replicates)
      each method
     
      output:
      val("${sc_dataset}-impact_technology")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_mixed_simulations.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{params.preProcess_dir}' '!{replicates}' '!{method}' '!{params.results_dir_impact_technology}' '!{params.ncores}' '!{baseDir}'
      ''' 
}

process ANALYSIS_BULK_MIRRORDB {

      label 'process_default'

      input:
      each sc_dataset
      each bulk_name 
      each method
     
      output:
      val("${sc_dataset}-impact_technology")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_mixed_simulations_real_dataset.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{bulk_name}' '!{params.data_dir_bulk}' '!{params.preProcess_dir}' '!{method}' '!{params.results_dir_impact_technology}' '!{params.ncores}' '!{baseDir}'
      ''' 
}

workflow impact_missing_cell_types {

  deconvolution = ANALYSIS_BULK_MISSING_CELL_TYPES(params.single_cell_list,
                                                   params.bulk_list,
                                                   params.method_list,
                                                   params.cell_types_to_exclude)

}

workflow simulation_impact_technology {

  simulations = SIMULATE_PSEUDOBULK_MIRRORDB(params.simulation_sc_dataset,
                                             params.simulation_n_cells,
						         params.simulation_n_samples,
						         params.sc_datasets_impact_technology, 
                                             params.replicates_simulation)

  deconvolution = ANALYSIS_PSEUDOBULK_MIRRORDB(params.simulation_sc_dataset,
                                               simulations, 
                                               params.method_list)

  deconvolution_real = ANALYSIS_BULK_MIRRORDB(params.simulation_sc_dataset,
                                              params.bulk_impact_technology,
                                              params.method_list)                                                                                      
}

workflow simulation_spillover {

  // spillover
  simulations = SIMULATE_BULK_SPILLOVER(params.simulation_sc_dataset,
                                        params.simulation_n_cells,
						    params.simulation_n_samples,
						    params.spillover_celltypes)
  
  deconvolution = ANALYSIS_SPILLOVER(simulations,
                                     params.spillover_celltypes,
                                     params.method_list,
                                     params.single_cell_list)

}

workflow simulation_unknown_content {
  // sensitivity to unknown content
  simulations = SIMULATE_BULK_UNKNOWN_CELL_TYPE(params.simulation_sc_dataset,
                                                params.simulation_n_cells,
							      params.simulation_n_samples,
							      params.fractions_unknown_cells,
                                                params.known_cell_types,
                                                params.unknown_cell_type, 
                                                params.replicates_unknown_content)
  
  deconvolution = ANALYSIS_BULK_UNKNOWN_CELL_TYPE(simulations,
                                             params.fractions_unknown_cells,
							   params.known_cell_types,
							   params.unknown_cell_type,
                                             params.method_list)
  
}

workflow simulation_resolution_analysis {
  // sensitivity to unknown content
  simulations = SIMULATE_BULK_RESOLUTION_ANALYSIS(params.simulation_sc_dataset,
                                                    params.simulation_n_cells,
							          params.simulation_n_samples,
							          params.cell_types_finer_res,
                                                    params.replicates_simulation)
  
  
  
  
  deconvolution = ANALYSIS_BULK_RESOLUTION_ANALYSIS(simulations,
                                                    params.cell_types_finer_res,
                                                    params.method_list)
  
}
