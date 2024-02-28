#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process PREPROCESS_SINGLE_CELL {

      input: 
      each sc_ds
      each method
      each ct_fractions
      each replicate

      output: 
      tuple val(sc_ds), 
            val(method),
            val(ct_fractions),
            val(replicate)

      shell:
      '''
      preprocessSingleCellNF.R '!{sc_ds}' '!{params.data_dir_sc}' '!{method}' '!{ct_fractions}' '!{replicate}' '!{params.preProcess_dir}'
      '''
}

process SIMULATE_BULK {
  
      //publishDir "${params.preProcess_dir}/pseudo_bulk", mode: 'copy'

      input:
      each simulation_n_cells
      each simulation_n_samples
      each simulation_scenario

      output:
      val("${params.simulation_sc_dataset}-ncells${simulation_n_cells}-nsamples${simulation_n_samples}-${simulation_scenario}")

      shell:
      '''
      /vol/omnideconv_input/benchmark/pipeline/bin/simulateBulkNF.R '!{params.simulation_sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{simulation_scenario}' '!{params.preProcess_dir}' '!{params.ncores}'
      '''
}

process ANALYSIS_BULK_MISSING_CELL_TYPES {

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
      analysisNF_missing_cell_type.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{bulk_ds}' '!{params.data_dir_bulk}' '!{method}' '!{cell_types_missing}' '!{params.results_dir_missing_cell_types}' '!{params.ncores}'
      ''' 
}

process SIMULATE_BULK_SPILLOVER {
      
      publishDir "${params.preProcess_dir}/pseudo_bulk_spillover", mode: 'copy'

      input:
      each sc_dataset
      each simulation_n_cells
      each simulation_n_samples
      val(cell_types)

      output:
      tuple val(sc_dataset),
            val("${sc_dataset}_spillover_sim"),
            val("${params.preProcess_dir}/pseudo_bulk_spillover/${sc_dataset}_spillover_sim")

      beforeScript 'chmod o+rw .'      
      shell:
      '''
      simulateBulkNF_spillover_analysis.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{cell_types}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_SPILLOVER {

      input:
      tuple val(sc_dataset), 
            val(sim_bulk_name), 
            val(sim_bulk_path)
      val cell_types
      each method 
     
      output:
      val("${sc_dataset}-spillover")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_spillover.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{method}' '!{cell_types}' '!{params.results_dir_spillover}' '!{params.ncores}'
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
      analysisNF_unknown_content.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{replicates}' '!{method}' '!{fraction_unknown_cell}' '!{cell_types}' '!{unknown_cell_type}' '!{params.results_dir_unknown_content}' '!{params.ncores}'
      ''' 
}

process SIMULATE_BULK_RESOLUTION_ANALYSIS {
  
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
            val("${params.preProcess_dir}/pseudo_bulk_resolution_newMethod/${sc_dataset}_resolution_analysis_sim"),
            val(replicates)
      
      beforeScript 'chmod o+rw .'     

      shell:
      '''
      simulateBulkNF_impact_cell_resolution.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{cell_types_fine}' '!{replicates}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_BULK_RESOLUTION_ANALYSIS {
  
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
      analysisNF_impact_cell_resolution.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{method}' '!{cell_types_fine}' '!{replicates}' '!{params.results_dir_resolution}' '!{params.ncores}'
      ''' 
}

process SIMULATE_PSEUDOBULK_MIRRORDB {
      
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
      simulateBulkNF_general.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{sc_dataset_list}' '!{replicates}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_PSEUDOBULK_MIRRORDB {

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
      analysisNF_mixed_simulations.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{sim_bulk_name}' '!{sim_bulk_path}' '!{params.preProcess_dir}' '!{replicates}' '!{method}' '!{params.results_dir_impact_technology}' '!{params.ncores}'
      ''' 
}

process ANALYSIS_BULK_MIRRORDB {

      input:
      each sc_dataset
      tuple val(bulk_name), 
            val(bulk_path)
      each method
     
      output:
      val("${sc_dataset}-impact_technology")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_mixed_simulations_real_dataset.R '!{sc_dataset}' '!{params.data_dir_sc}' '!{bulk_name}' '!{bulk_path}' '!{params.preProcess_dir}' '!{method}' '!{params.results_dir_impact_technology}' '!{params.ncores}'
      ''' 
}

process CREATE_SIGNATURE {

      input:
      each sc_ds
      each bulk_ds
      each method 
      each ct_fractions
      each replicate
      val run_preprocessing

      output:
      tuple val(sc_ds), 
            val(bulk_ds),
            val(method),
            val(replicate), 
            val(ct_fractions)

      beforeScript 'chmod o+rw .'

      shell:
      '''
      computeSignaturesNF.R '!{sc_ds}' '!{params.data_dir_sc}' '!{bulk_ds}' '!{params.data_dir_bulk}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
      ''' 
}

process CREATE_SIGNATURE_PREPROCESSED {

      input:
      tuple val(sc_ds), 
            val(method),
            val(ct_fractions),
            val(replicate)
      val sc_path
      each bulk_ds
      val run_preprocessing

      output:
      tuple val(sc_ds), 
            val(bulk_ds),
            val(method),
            val(replicate), 
            val(ct_fractions)

      beforeScript 'chmod o+rw .'

      shell:
      '''
      computeSignaturesNF.R '!{sc_ds}' '!{sc_path}' '!{bulk_ds}' '!{params.data_dir_bulk}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
      ''' 
}

process CREATE_SIGNATURE_RECTANGLE {

      input:
      tuple path(sc_matrix), 
            path(sc_anno), 
            path(sc_batch), 
            val(sc_ds), 
            val(sc_norm), 
            val(replicate), 
            val(ct_fractions)
      val bulk_dir
      each bulk_ds
      each bulk_norm
      each method 
      val run_preprocessing

      output:
      tuple path(sc_matrix), 
            path(sc_anno), 
            path(sc_batch), 
            val(sc_ds), 
            val(sc_norm),
            val(bulk_ds),
            val(bulk_norm),
            val(replicate), 
            val(ct_fractions),
            val(method)

      beforeScript 'chmod o+rw .'

      shell:
      '''
      /vol/omnideconv_input/benchmark/pipeline/bin/rectangle/computeSignaturesNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
      ''' 
}

process CREATE_SIGNATURE_FOR_SIMULATION {

      input:
      tuple path(sc_matrix), 
            path(sc_anno), 
            path(sc_batch), 
            val(sc_ds), 
            val(sc_norm), 
            val(replicate), 
            val(ct_fractions)
      val bulk_dir
      each bulk_ds
      each bulk_norm
      val cell_types
      each method 
      val run_preprocessing


      output:
      tuple path(sc_matrix), 
            path(sc_anno), 
            path(sc_batch), 
            val(sc_ds), 
            val(sc_norm),
            val(bulk_ds),
            val(bulk_norm),
            val(replicate), 
            val(ct_fractions),
            val(method)

      beforeScript 'chmod o+rw .'

      shell:
      '''
      /vol/omnideconv_input/benchmark/pipeline/bin/computeSignaturesNF_simulation.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{cell_types}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
      ''' 
}

process DECONVOLUTE { 

	input:
      tuple val(sc_ds), 
            val(bulk_ds),
            val(method),
            val(replicate), 
            val(ct_fractions)
      val sc_path
	val run_preprocessing
  
	output:
	tuple val(method),
	      val(sc_ds), 
	      val(bulk_ds), 
	      val(replicate),
	      val(ct_fractions)
  
	shell:
	'''
	runDeconvolutionNF.R '!{sc_ds}' '!{sc_path}' '!{bulk_ds}' '!{params.data_dir_bulk}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.species_sc}' '!{params.ncores}'
	''' 
}

process COMPUTE_METRICS { 

	input:
	tuple val(method),
	      val(sc_ds), 
	      val(bulk_ds), 
	      val(replicate),
	      val(ct_fractions)
  
	output:
	tuple val(method),
	      val(sc_ds), 
	      val(bulk_ds), 
	      val(replicate),
	      val(ct_fractions)
  
	shell:
	'''
	computeMetricsNF.R '!{sc_ds}' '!{params.data_dir_sc}' '!{bulk_ds}' '!{params.data_dir_bulk}' '!{method}' '!{replicate}' '!{ct_fractions}' '!{params.results_dir_general}'
	''' 
}

workflow simulation {
  
  // Normal simulation
  
  simulations = SIMULATE_BULK(params.simulation_n_cells,
							                params.simulation_n_samples,
							                params.simulation_scenario
  )
  
  sc_files = Channel.fromList(create_file_list_sc(params.data_dir_sc, 
                                                  params.single_cell_list, 
                                                  params.single_cell_norm, 
                                                  'false'))
  
  
  signature = CREATE_SIGNATURE(sc_files,
                               "${params.preProcess_dir}/pseudo_bulk",
                               simulations.collect(),
                               params.simulation_pseudobulk_norm,
                               params.simulation_cell_types,
                               params.method_list,
                               'false')
  
  deconvolution = DECONVOLUTE(signature, 
                             "${params.preProcess_dir}/pseudo_bulk",
                             'false')
  
  metrics = COMPUTE_METRICS(deconvolution, "${params.preProcess_dir}/pseudo_bulk")
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
						         params.datasets_impact_technology, 
                                             params.replicates_simulation)

  //deconvolution = ANALYSIS_BULK_MIRRORDB(params.simulation_sc_dataset,
  //                                       ['vanderbilt_lung', '/nfs/data/omnideconv_benchmarking_clean/data/Tumor'], 
  //                                       params.method_list)          

  deconvolution = ANALYSIS_PSEUDOBULK_MIRRORDB(params.simulation_sc_dataset,
                                               simulations, 
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
                                     params.method_list)

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

workflow subsampling {

  replicates = Channel.of(1..params.replicates).collect()
    
  preprocess = PREPROCESS_SINGLE_CELL(params.single_cell_list,
                                      params.method_list,
                                      params.ct_fractions, 
                                      replicates
  )
  
  signature = CREATE_SIGNATURE_PREPROCESSED(preprocess,
                                            params.preProcess_dir,
                                            params.bulk_list,
                                            'true'
  )  
  
  deconvolution = DECONVOLUTE(signature,
                              params.preProcess_dir,
                              'true')

  metrics = COMPUTE_METRICS(deconvolution)
  
}

workflow {

  signature = CREATE_SIGNATURE(params.single_cell_list,
                               params.bulk_list,
                               params.method_list,
                               0,
                               0,
                               'false'
  )  
  
  deconvolution = DECONVOLUTE(signature,
                              params.data_dir_sc,      
                              'false'
  )

  metrics = COMPUTE_METRICS(deconvolution) 
}