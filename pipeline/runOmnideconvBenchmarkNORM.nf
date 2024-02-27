#!/usr/bin/env nextflow
nextflow.enable.dsl=2

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

def create_file_list_sc(basepath, ds_list, norm_list, do_preprocessing) {
    def file_list = []
    for (sc_ds in ds_list) {
        for (sc_norm in norm_list) {
            def anno = "${basepath}/${sc_ds}/celltype_annotations.rds"
            def batch = "${basepath}/${sc_ds}/batch.rds"
            def matrix = "${basepath}/${sc_ds}/matrix_"
            
            if (sc_norm == 'counts') {
                matrix += "counts.rds"
            } else {
                matrix += "norm_counts.rds"
            }
            
            if (do_preprocessing == 'true'){
                file_list << [matrix, anno, batch, sc_ds, sc_norm]
            }else{
                file_list << [matrix, anno, batch, sc_ds, sc_norm, 0, 0]
            }
        }
    }
    return file_list
}


process PREPROCESS_SINGLE_CELL {

      input: 
      tuple path(sc_matrix), 
      path(sc_anno), 
      path(sc_batch), 
      val(sc_ds), 
      val(sc_norm)
      each ct_fractions
      each replicate
      path preProcess_dir

      output: 
      tuple path("${preProcess_dir}/${sc_ds}_${sc_norm}_perc${ct_fractions}_rep${replicate}/matrix_subsampled.rds"), 
      path("${preProcess_dir}/${sc_ds}_${sc_norm}_perc${ct_fractions}_rep${replicate}/celltype_annotations.rds"), 
      path("${preProcess_dir}/${sc_ds}_${sc_norm}_perc${ct_fractions}_rep${replicate}/batch.rds"), 
      val(sc_ds), 
      val(sc_norm), 
      val(replicate),
      val(ct_fractions)

      shell:
      '''
      preprocessSingleCellNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{ct_fractions}' '!{replicate}' '!{preProcess_dir}'
      '''
}

process SIMULATE_BULK {
  
      publishDir "${params.preProcess_dir}/pseudo_bulk", mode: 'copy'

      input:
      each simulation_n_cells
      each simulation_n_samples
      each simulation_scenario

      output:
      val("${params.simulation_sc_dataset}-ncells${simulation_n_cells}-nsamples${simulation_n_samples}-${simulation_scenario}")

      shell:
      '''
      simulateBulkNF.R '!{params.simulation_sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{simulation_scenario}' '!{params.preProcess_dir}' '!{params.ncores}'
      '''
}

process SIMULATE_BULK_SPILLOVER {
  
      publishDir "${params.preProcess_dir}/pseudo_bulk_spillover", mode: 'copy'

      input:
      val simulation_n_cells
      val spillover_samples_per_cell
      val cell_types

      output:
      val("${params.simulation_sc_dataset}-spillover")

      shell:
      '''
      simulateBulkNF_spillover_analysis.R '!{params.simulation_sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{spillover_samples_per_cell}' '!{params.cell_types}' '!{params.preProcess_dir}' '!{params.ncores}' 
      '''
}

process ANALYSIS_SPILLOVER {
  
      publishDir "${params.preProcess_dir}/pseudo_bulk_resolution", mode: 'copy'

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
      val("${params.simulation_sc_dataset}-spillover")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_spillover.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{cell_types}' '!{params.results_dir_spillover}' '!{run_preprocessing}' '!{params.ncores}'
      ''' 
}


process SIMULATE_BULK_SENSITIVITY {
  
      publishDir "${params.preProcess_dir}/pseudo_bulk_sensitivity", mode: 'copy'

      input:
      each simulation_n_cells
      val simulation_n_samples                    
      each fraction_unknown_cell
      val known_cell_types
      val unknown_cell_type

      output:
      val("${params.simulation_sc_dataset}-fraction_${unknown_cell_type}_${fraction_unknown_cell}-unknown")

      shell:
      '''
      simulateBulkNF_unknown_content_analysis.R '!{params.simulation_sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{params.fraction_unknown_cell}' '!{params.preProcess_dir}' '!{params.ncores}' '!{params.known_cell_types}' '!{params.unknown_cell_type}' 
      '''
}


process ANALYSIS_SENSITIVITY {
  
      publishDir "${params.preProcess_dir}/pseudo_bulk_resolution", mode: 'copy'

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
      val("${params.simulation_sc_dataset}-fraction_${unknown_cell_type}_${fraction_unknown_cell}-unknown")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_unknown_content.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{cell_types}' '!{params.results_dir_unknown_content}' '!{run_preprocessing}' '!{params.ncores}'
      ''' 
}


process SIMULATE_BULK_RESOLUTION_ANALYSIS {
  
      publishDir "${params.preProcess_dir}/pseudo_bulk_resolution", mode: 'copy'

      input:
      each simulation_n_cells
      val simulation_n_samples                   
      each cell_types_finer_res
      output:
      val("${params.simulation_sc_dataset}-resolution_analysis")

      shell:
      '''
      simulateBulkNF_impact_cell_resolution.R '!{params.simulation_sc_dataset}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{params.fraction_unknown_cell}' '!{params.preProcess_dir}' '!{params.ncores}' '!{params.cell_types_finer_res}' 
      '''
}

process ANALYSIS_RESOLUTION {
  
      publishDir "${params.preProcess_dir}/pseudo_bulk_resolution", mode: 'copy'

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
      val("${params.simulation_sc_dataset}-resolution_analysis")
      
      beforeScript 'chmod o+rw .'
      
      shell:
      '''
      analysisNF_impact_cell_resolution.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{cell_types}' '!{params.results_dir_resolution}' '!{run_preprocessing}' '!{params.ncores}'
      ''' 
}

process CREATE_SIGNATURE {

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
      computeSignaturesNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
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
      computeSignaturesNF_simulation.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{cell_types}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
      ''' 
}

process DECONVOLUTE { 

	input:
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
	val bulk_dir
	val run_preprocessing
  
	output:
	tuple val(method),
	      val(sc_ds), 
	      val(sc_norm), 
	      val(bulk_ds), 
	      val(bulk_norm), 
	      val(replicate),
	      val(ct_fractions)
  
	shell:
	'''
	runDeconvolutionNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.species_sc}' '!{params.ncores}' 
	''' 
}

process COMPUTE_METRICS { 

	input:
	tuple val(method),
	      val(sc_ds), 
	      val(sc_norm), 
	      val(bulk_ds), 
	      val(bulk_norm), 
	      val(replicate),
	      val(ct_fractions)
	val bulk_dir
  
	output:
	tuple val(method),
	      val(sc_ds), 
	      val(sc_norm), 
	      val(bulk_ds), 
	      val(bulk_norm), 
	      val(replicate),
	      val(ct_fractions)
  
	shell:
	'''
	computeMetricsNF.R '!{sc_ds}' '!{sc_norm}' '!{bulk_dir}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{replicate}' '!{ct_fractions}' '!{params.results_dir_general}'
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

workflow other_simulation {
  // spillover
  
  simulations = SIMULATE_BULK_SPILLOVER(params.simulation_n_cells,
		                            params.spillover_samples_per_cell,
						    params.spillover_celltypes
  )
  
  sc_files = Channel.fromList(create_file_list_sc(params.data_dir_sc, 
                                                  params.single_cell_list, 
                                                  params.single_cell_norm, 
                                                  'false'))
  
  
  deconvolution = ANALYSIS_SPILLOVER(sc_files,
                                     "${params.preProcess_dir}/pseudo_bulk_spillover",
                                     simulations.collect(),
                                     params.simulation_pseudobulk_norm,
                                     params.spillover_celltypes,
                                     params.method_list,
                                     'false')
  
  
  // sensitivity to unknown content
  
  simulations = SIMULATE_BULK_SENSITIVITY(params.simulation_n_cells,
						      params.simulation_n_samples,
							params.fraction_unknown_cell,
							params.known_cell_types,
							params.unknown_cell_type
  )
  
  sc_files = Channel.fromList(create_file_list_sc(params.data_dir_sc, 
                                                  params.single_cell_list, 
                                                  params.single_cell_norm, 
                                                  'false'))
  
  
  deconvolution = ANALYSIS_SENSITIVITY(sc_files,
                                       "${params.preProcess_dir}/pseudo_bulk_sensitivity",
                                       simulations.collect(),   // I think this will need ot be edited but not 100% sure how 
                                       params.simulation_pseudobulk_norm,
                                       params.known_cell_types,
                                       params.method_list,
                                       'false')
  

  // resolution analysis
  
  simulations = SIMULATE_BULK_RESOLUTION_ANALYSIS(params.simulation_n_cells,
							        params.simulation_n_samples,
							        params.cell_types_finer_res
  )
  
  sc_files = Channel.fromList(create_file_list_sc(params.data_dir_sc, 
                                                  params.single_cell_list, 
                                                  params.single_cell_norm, 
                                                  'false'))
  
  
  deconvolution = ANALYSIS_RESOLUTION(sc_files,
                                      "${params.preProcess_dir}/pseudo_bulk_resolution",
                                      simulations.collect(),
                                      params.simulation_pseudobulk_norm,
                                      params.known_cell_types,
                                      params.method_list,
                                      'false')
  
  
}

workflow subsampling {
  sc_files = Channel.fromList(create_file_list_sc(params.data_dir_sc, 
                                                  params.single_cell_list, 
                                                  params.single_cell_norm, 
                                                  'true'))

  replicates = Channel.of(1..params.replicates).collect()
    
  preprocess = PREPROCESS_SINGLE_CELL(sc_files,
                                      params.ct_fractions, 
                                      replicates, 
                                      params.preProcess_dir,
                                      'true'
  )
  
  signature = CREATE_SIGNATURE(preprocess,
                               params.data_dir_bulk,
                               params.bulk_list,
                               params.bulk_norm,
                               params.method_list,
                               'true'
  )  
  
  deconvolution = DECONVOLUTE(signature, 
                              params.data_dir_bulk,
                              'true')

  metrics = COMPUTE_METRICS(deconvolution, params.data_dir_bulk)
  
}

workflow {
  sc_files = Channel.fromList(create_file_list_sc(params.data_dir_sc, 
                                                  params.single_cell_list, 
                                                  params.single_cell_norm, 
                                                  'false'))
  preprocess = Channel.empty()
  signature = CREATE_SIGNATURE(sc_files,
                               params.data_dir_bulk,
                               params.bulk_list,
                               params.bulk_norm,
                               params.method_list,
                               'false'
  )  
  
  deconvolution = DECONVOLUTE(signature, params.data_dir_bulk, 'false')

  metrics = COMPUTE_METRICS(deconvolution, params.data_dir_bulk)
}