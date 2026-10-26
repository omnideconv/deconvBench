# Changelog

## 2.0.0

The pipeline is rebuilt as a modular, file-based Nextflow pipeline (Nextflow 25.10). Results of the existing methods are unchanged except for the fixes listed below; the published result layout is kept, so the scripts in `visualisation/` work as before.

**Pipeline**
- One process per module, workflows in `main.nf`; datasets are process inputs, so `-resume` re-runs what depends on a changed file.
- Signatures that depend only on the reference are built once per reference; pseudo-bulks are simulated once and shared by all methods; all simulations are seeded.
- One cluster config (`conf/norm.config`) with a profile per image and memory per method, short run configs in `configs/`, launchers `scripts/run.sh` and `scripts/run_profiles.sh`.
- The main workflow deconvolves bulk datasets without ground truth and skips their metrics.

**New**
- Methods: rectangle (rectanglepy 1.7.0), the hierarchical MuSiC (`music-hc`) and SECRET.
- Split of the single-cell reference between simulation and deconvolution, by donor or within every cell type (`split_reference`, `split_by`).
- Continuous integration on GitHub: all workflows and every method in its image on small test data (`data/ci/`).
- Images as `.sif` files in the GitHub container registry.

**Fixes that change results**
- Unknown content: the SimBu whitelist was empty, so pseudo-bulks contained every cell type of the dataset.
- Missing cell types: Scaden's cell-type names came back with dots, so metrics covered only one-word names.
- AutoGeneS: `+`, `-`, `(` and `)` in cell-type names were not translated back.
- Resolution analysis: a normal or coarse cell type with a single fine cell type of another name kept the fine name in the ground truth.
- rectangle: the number of cores was not passed through.
- Seeds for simulations that had none.

## 1.1.0

Rectangle benchmark.

## 1.0.0

First release.
