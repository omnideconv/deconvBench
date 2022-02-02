#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process setup {

	publishDir params.results_dir, mode: params.publishDirMode

	output:
	path('*.tsv'), emit: runtime

	"""
	printf "method\tsc_dataset\trnaseq_dataset\tstep\truntime_user\truntime_system\truntime_elapsed\n" > runtimes.tsv
	""" 
}

process createSignature {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each sc
	each rnaseq
	each mode 
	path rnaseq_data
	path sc_data
	path mapping_sheet
	path runtime

	output:
	tuple val(sc), val(rnaseq), val(mode), path('*.rds'), emit: signature
	path(runtime), emit: runtime_sheet

  shell:
	'''
	#echo '$sc $sc_data $rnaseq_data $mode $rnaseq'
	#computeSignaturesNF.R '$sc_data' '$sc' '$rnaseq_data' '$rnaseq' '$mode' '$mapping_sheet'
	
	TIME="$(computeSignaturesNF.R '!{sc_data}' '!{sc}' '!{rnaseq_data}' '!{rnaseq}' '!{mode}' '!{mapping_sheet}')"
	FINAL="$(echo $TIME | sed -e 's/^.*user system elapsed //p' | tail -n 1)"
	echo !{mode}'\t'!{sc}'\t'!{rnaseq}'\tbuild_model\t'$FINAL >> '!{runtime}'
	''' 
}

process deconvolute { 

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	tuple val(sc), val(rnaseq), val(mode), val(signature)
	path rnaseq_data
	path sc_data
	path mapping_sheet
	path runtime

	output:
	path('*.rds'), emit: deconvolution
	path(runtime), emit: runtime

	shell:
	'''
	#echo '$sc $rnaseq $signature $sc_data $rnaseq_data'
	#runDeconvolutionNF.R '$sc_data' '$sc' '$rnaseq_data' '$rnaseq' '$mode' '$signature' '$mapping_sheet'
	TIME="$(runDeconvolutionNF.R '!{sc_data}' '!{sc}' '!{rnaseq_data}' '!{rnaseq}' '!{mode}' '!{signature}' '!{mapping_sheet}')"
	FINAL="$(echo $TIME | sed -e 's/^.*user system elapsed //p' | tail -n 1)"
	echo !{mode}'\t'!{sc}'\t'!{rnaseq}'\tdeconvolute\t'$FINAL >> '!{runtime}'
	''' 
}

process benchmarkingPlots {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	val deconvolution
	path facs
	val dataset_name
	path mapping_sheet

	
	output:
	//stdout emit: plot
	path '*.jpeg', emit: plotjpeg
	path '*.pdf', emit: plotpdf

	"""
	echo '$deconvolution'
	plotGtruthHoekNF.R '$deconvolution' '$facs' '$dataset_name' '$mapping_sheet'
	""" 
}

process simulateBulk {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each sc
	each scenario 
	path sc_data
	path mapping_sheet
	val celltypes

	output:
	path('*.rds'), emit: bulk
	//path('*.pdf'), emit: correlation

	"""
	#echo '$sc $sc_data $scenario'
	simulateBulkNF.R '$sc_data' '$sc' '$mapping_sheet' '$scenario' 'blood' 'Homo sapiens' '$celltypes'
	""" 
}

process deconvoluteSimulatedBulk {

	publishDir params.results_dir, mode: params.publishDirMode

	input:
	each sc
	each method
	each bulk
	path sc_data
	path mapping_sheet

	output:
	path('*.rds'), emit: deconvolution
	path('*.pdf'), emit: correlation

	"""
	#echo '$sc $sc_data $method $bulk'
	deconvoluteSimBulkNF.R '$sc_data' '$sc' '$bulk' '$method' '$mapping_sheet'
	""" 
}


workflow{
  setup()
  runtime = setup.out.runtime //is there an easier way?
  single_cell = Channel.fromList(params.single_cell_list)
  rna_seq = Channel.fromList(params.rnaseq_list)
  methods = params.method_list
  scenarios = params.simulation_scenarios
  rnaseq_data = Channel.fromPath(params.data_dir_rnaseq).collect()
  sc_data = Channel.fromPath(params.data_dir_sc).collect()
  remapping = Channel.fromPath(params.mapping_sheet).collect()
  //createSignature(single_cell, rna_seq, methods, rnaseq_data, sc_data, remapping, runtime)
  //deconvolute(createSignature.out.signature, rnaseq_data, sc_data, remapping, runtime)
  //deconv = deconvolute.out.deconvolution
  //deconvolute.out.deconvolution.view()
  //deconvolute.out.runtime.view()
  //hoek_samples = deconv.filter{ ds -> ds =~/_hoek/ }                      das
  //hoek_samples_list = hoek_samples.toList()                               hier
  //benchmarkingPlots(hoek_samples_list, rnaseq_data, rna_seq, remapping)   nicht
  //benchmarkingPlots(deconv.toList(), rnaseq_data, rna_seq, remapping)
  //benchmarkingPlots.out.plot.view()
  celltypes = Channel.fromList(["B cell", "T cell CD4+", "T cell CD8+"])
  simulateBulk(single_cell, scenarios, sc_data, remapping, celltypes)
  sims = simulateBulk.out.bulk.toList()
  deconvoluteSimulatedBulk(single_cell, methods, sims, sc_data, remapping)
  deconvSim = deconvoluteSimulatedBulk.out.deconvolution
  spillover_samples = deconvSim.filter{ str -> str =~/_spillover/ }.toList()
}