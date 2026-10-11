include { splitSetting; splitFolder } from '../../functions'

// Split a single-cell dataset by donor into a part for the simulation of pseudo-bulks and a part
// that the methods use as reference, so that no donor contributes to both.
process SPLIT_REFERENCE {

      tag "${sc_name}"
      label 'process_default'

      // the tables that describe the split; the two parts stay in the work directory
      publishDir path: { "${params.preProcess_dir}/reference_split/${splitFolder(sc_name)}" }, mode: 'copy', pattern: 'split/*', saveAs: { filename -> filename.replaceFirst('^split/', '') }

      input:
      tuple val(sc_name), path(sc_dir, stageAs: 'sc/*')

      output:
      tuple val(sc_name), path("simulation/${sc_name}"), emit: simulation
      tuple val(sc_name), path("deconvolution/${sc_name}"), emit: deconvolution
      path 'split/*', emit: split

      script:
      """
      splitReferenceNF.R '${sc_name}' 'sc' '${splitSetting('split_fraction_simulation')}' '${splitSetting('split_min_cells')}' '${splitSetting('split_min_cells_action')}' '${splitSetting('split_seed')}' '${splitSetting('split_reference_cells_fraction')}' '.'
      """
}
