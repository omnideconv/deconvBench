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
            
            if (do_preprocessing){
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
    /vol/omnideconv_input/benchmark/pipeline/bin/preprocessSingleCellNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{ct_fractions}' '!{replicate}' '!{preProcess_dir}'
    '''
}

process SIMULATE_BULK {
  
  input:
  each simulation_n_cells
  each simulation_n_samples
  each simulation_scenario
  
  output:
  tuple val("${simulation_sc_dataset}_ncells${simulation_n_cells}_nsamples${simulation_n_samples}_${simulation_scenario}")
        val("${simulation_sc_norm}")
  //path("${preProcess_dir}/pseudo_bulk/${simulation_sc_dataset}_${simulation_sc_norm}_ncells${simulation_n_cells}_nsamples${simulation_n_samples}_${simulation_scenario}/pseudobulk.rds")
  //path("${preProcess_dir}/pseudo_bulk/${simulation_sc_dataset}_${simulation_sc_norm}_ncells${simulation_n_cells}_nsamples${simulation_n_samples}_${simulation_scenario}/true_fractions.rds")
        

  shell:
  '''
  /vol/omnideconv_input/benchmark/pipeline/bin/simulateBulkNF.R '!{params.simulation_sc_dataset}' '!{params.simulation_sc_norm}' '!{params.data_dir_sc}' '!{simulation_n_cells}' '!{simulation_n_samples}' '!{simulation_scenario}' '!{params.preProcess_dir}' '!{params.ncores}'
  '''
}

process CREATE_SIGNATURE_FROM_SIMULATION {

  input:
  tuple path(sc_matrix), 
	      path(sc_anno), 
	      path(sc_batch), 
	      val(sc_ds), 
	      val(sc_norm), 
	      val(replicate), 
	      val(ct_fractions)
  each pseudobulk
  each pseudobulk_norm
  each method

  output:
  tuple path(sc_matrix), 
        path(sc_anno), 
        path(sc_batch), 
        val(sc_ds), 
        val(sc_norm),
        val(pseudobulk),
        val(pseudobulk_norm),
        val(replicate), 
        val(ct_fractions),
        val(method)

  beforeScript 'chmod o+rw .'

	shell:
	'''
	/vol/omnideconv_input/benchmark/pipeline/bin/computeSignaturesNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{params.preProcess_dir}/pseudo_bulk' '!{pseudobulk}' '!{pseudobulk_norm}' '!{method}' '!{params.results_dir_general}' 'false' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
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
	each bulk_ds
	each bulk_norm
	each method 
	val(run_preprocessing)
  
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
	/vol/omnideconv_input/benchmark/pipeline/bin/computeSignaturesNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{params.data_dir_bulk}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{ct_fractions}' '!{params.ncores}'
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
	val(run_preprocessing)
  
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
	/vol/omnideconv_input/benchmark/pipeline/bin/runDeconvolutionNF.R '!{sc_matrix}' '!{sc_anno}' '!{sc_batch}' '!{sc_ds}' '!{sc_norm}' '!{params.data_dir_bulk}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{params.results_dir_general}' '!{run_preprocessing}' '!{replicate}' '!{params.species_sc}' '!{ct_fractions}' '!{params.ncores}' 
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
	/vol/omnideconv_input/benchmark/pipeline/bin/computeMetricsNF.R '!{sc_ds}' '!{sc_norm}' '!{params.data_dir_bulk}' '!{bulk_ds}' '!{bulk_norm}' '!{method}' '!{replicate}' '!{ct_fractions}' '!{params.results_dir_general}'
	''' 
}


workflow simulation {
  
  simulations = SIMULATE_BULK(params.simulation_n_cells,
							                params.simulation_n_samples,
							                params.simulation_scenario
  )


  
  //signature = CREATE_SIGNATURE_FROM_SIMULATION()
  
  //deconvolution = DECONVOLUTE(signature)
  
  //metrics = COMPUTE_METRICS(deconvolution)
  
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
                                      params.preProcess_dir
  )
  
  signature = CREATE_SIGNATURE(preprocess,
                               params.bulk_list,
                               params.bulk_norm,
                               params.method_list,
                               'true'
  )  
  
  deconvolution = DECONVOLUTE(signature, 'true')

  metrics = COMPUTE_METRICS(deconvolution)
  
}

workflow {
  sc_files = Channel.fromList(create_file_list_sc(params.data_dir_sc, 
                                                  params.single_cell_list, 
                                                  params.single_cell_norm, 
                                                  'false'))
  preprocess = Channel.empty()
  signature = CREATE_SIGNATURE(preprocess,
                               params.bulk_list,
                               params.bulk_norm,
                               params.method_list,
                               'false'
  )  
  
  deconvolution = DECONVOLUTE(signature, 'false')

  metrics = COMPUTE_METRICS(deconvolution)
}
