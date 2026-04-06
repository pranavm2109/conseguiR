# conseguiR External API Overview

This file describes the intended user-callable external API in plain English.

The external API is designed as a sequence of bundle-returning functions:

1. raw inputs -> score bundles
2. score bundles -> scored graph bundle
3. scored graph bundle -> diffusion bundle
4. diffusion bundle -> selected subgraph bundle
5. selected subgraph bundle -> plot bundle

Every external function returns a bundle with:

- `bundle_type`: bundle label
- `objects`: in-memory tables and graphs
- `output_paths`: saved file paths
- `config`: arguments used to build the result

Many bundles also expose their main objects at top level for convenience.

## Shared Terminology

Across the external API docs, the following terms are used consistently:

- `bundle`: a named result object with `objects`, `output_paths`, and `config`
- `wrapper-level input`: a parameter exposed directly by the user-facing R
  function
- `passthrough argument`: a parameter forwarded through a list such as
  `extra_args`, `dndscv_args`, or `fishhook_args`
- `output path`: a file path recorded in the returned bundle for reuse by a
  downstream stage
- `stage`: one computational step in the end-to-end pipeline

## Validation

### `validate_inputs()`

Purpose:
- validate raw user inputs before scoring

Main inputs:
- `gwas_sumstats`
- `somatic_maf`
- `reg_ref_path`
- `epigenomic_tracks` or `epigenomic_track_dir`

Main outputs:
- validated GWAS table
- validated somatic MAF table
- validated regulatory GRanges
- validated epigenomic input summary

Returns:
- validation bundle

## Germline Scoring

For a more detailed MAGMA option map, see:

- `scripts/Externals/R/MAGMA_ARGUMENT_COVERAGE.md`

### `run_germline_gene_scoring()`

Purpose:
- run MAGMA-based germline scoring for genes

Main inputs:
- `gwas_sumstats`
- `gene_loc_path`
- `reference_bfile`
- separate `step1_args` and `step2_args`

Stage split:
- MAGMA step 1 = annotation
- MAGMA step 2 = gene analysis

Typical `step1_args` entries:
- `annotation_window`
- `filter_path`
- `ignore_strand`
- `nonhuman`
- `extra_args`

Typical `step2_args` entries:
- `gene_model`
- `genes_only`
- `pval_use`
- `pval_duplicate`
- `bfile_synonyms`
- `bfile_synonym_dup`
- `extra_args`

Core wrapper-level MAGMA step 2 inputs:
- `sample_size`
- `sample_size_col`

Main outputs:
- gene germline score table
- internal MAGMA pipeline result

Returns:
- germline gene score bundle

### `run_germline_regulatory_scoring()`

Purpose:
- run MAGMA-based germline scoring for regulatory elements

Main inputs:
- `gwas_sumstats`
- `reg_loc_path`
- `reference_bfile`
- separate `step1_args` and `step2_args`

Stage split:
- MAGMA step 1 = annotation
- MAGMA step 2 = gene analysis

The same MAGMA stage-specific structure is used here, but for the regulatory
run rather than the gene run.

Main outputs:
- regulatory germline score table
- internal MAGMA pipeline result

Returns:
- germline regulatory score bundle

### `prepare_germline_scores()`

Purpose:
- orchestrate both MAGMA germline runs

Main inputs:
- `gwas_sumstats`
- `reference_bfile`
- `gene_loc_path`
- `reg_loc_path`
- separate parameter lists for gene step 1, gene step 2, reg step 1, reg step 2

Concretely, this function exposes MAGMA stage customization twice:
- `gene_step1_args`
- `gene_step2_args`
- `reg_step1_args`
- `reg_step2_args`

So the user can tune the gene MAGMA run and the regulatory MAGMA run
independently.

Main outputs:
- `gene_scores`
- `reg_scores`
- underlying gene and reg result bundles

Returns:
- germline score bundle

## Somatic Scoring

For more detailed argument maps, see:

- `scripts/Externals/R/DNDSCV_ARGUMENT_COVERAGE.md`
- `scripts/Externals/R/FISHHOOK_ARGUMENT_COVERAGE.md`

### `run_somatic_gene_scoring()`

Purpose:
- run dndscv-based somatic scoring for genes

Main inputs:
- `maf`
- `refdb`
- `cv`
- `max_muts_per_gene_per_sample`
- `max_coding_muts_per_sample`
- `dndscv_args`

Main outputs:
- somatic gene score table

Returns:
- somatic gene score bundle

### `run_somatic_regulatory_scoring()`

Purpose:
- run fishHook-based somatic scoring for regulatory elements

Main inputs:
- `maf`
- `reg_ref_path`
- `eligible_gr`
- `fishhook_covariates`
- `fishhook_covariate_data`
- `idcol`
- `fishhook_args`

Main outputs:
- somatic regulatory score table

Returns:
- somatic regulatory score bundle

### `prepare_somatic_scores()`

Purpose:
- orchestrate both somatic scoring paths

Main inputs:
- `maf`
- `refdb`
- `reg_ref_path`
- wrapper-level dndscv settings
- wrapper-level fishHook settings

Main outputs:
- `gene_scores`
- `reg_scores`
- underlying gene and reg result bundles

Returns:
- somatic score bundle

## Epigenomic Scoring

### `prepare_epigenomic_scores()`

Purpose:
- compute regulatory-element epigenomic z scores from bigWig signal tracks

Main inputs:
- `reg_ref_path`
- `bw_files` or `track_dir`
- `exclude_patterns`
- `min_tracks`
- `transform`

Main outputs:
- regulatory epigenomic score table
- optional diagnostics

Returns:
- epigenomic score bundle

## Scored Graph Construction

### `build_scored_gene_reg_graph()`

Purpose:
- impose all score modalities onto the backend no-score gene-reg graph

Main inputs:
- no-score graph path or graph object
- germline score bundle or explicit gene/reg germline score tables
- somatic score bundle or explicit gene/reg somatic score tables
- epigenomic score bundle or explicit regulatory epigenomic score table

Behavior:
- genes receive somatic and germline scores
- genes get `epigenomic_score = 0`
- regulatory elements receive somatic, germline, and epigenomic scores
- missing scores default to zero

Main outputs:
- scored `igraph`
- scored node table
- scored edge table

Returns:
- scored gene-reg graph bundle

## Diffusion

For a more detailed argument map, see:

- `scripts/Externals/R/DIFFUSION_ARGUMENT_COVERAGE.md`

### `run_gene_reg_diffusion()`

Purpose:
- run the Python-backed diffusion step on the scored gene-reg graph

Main inputs:
- scored graph bundle or explicit node/edge paths
- diffusion hyperparameters

Main outputs:
- full diffusion table
- top genes diffusion table

Returns:
- diffusion bundle

## Subgraph Calling

For a more detailed argument map, see:

- `scripts/Externals/R/SUBGRAPH_ARGUMENT_COVERAGE.md`

### `call_selected_subgraph()`

Purpose:
- run the Python-backed cardinality-constrained subgraph selection step

Main inputs:
- diffusion bundle or diffusion path
- gene-gene node and edge paths
- target subgraph size
- solver settings

Main outputs:
- selected subgraph node table
- selected subgraph edge table
- selected subgraph summary table
- GraphML path

Returns:
- selected subgraph bundle

## Plotting

For a more detailed argument map, see:

- `scripts/Externals/R/PLOTTING_ARGUMENT_COVERAGE.md`

### `plot_scores()`

Purpose:
- create a rank plot for one-tailed outputs and a volcano plot for two-tailed outputs

Main inputs:
- score or diffusion bundle, or explicit table
- bundle component selector such as `gene_scores`, `reg_scores`, `all_genes`, or `top_genes`
- tail direction (`one_tailed` or `two_tailed`)
- optional feature labels to highlight

Main outputs:
- ggplot object
- optional saved plot file

Returns:
- plot bundle

### `plot_selected_subgraph()`

Purpose:
- build a visualization bundle and optionally save a viewable figure

Main inputs:
- selected subgraph bundle or explicit node/edge/summary inputs
- bundle save prefix
- figure `file_path`
- plot title and layout settings

Main outputs:
- ggplot object
- visualization bundle
- optional saved plot and saved bundle tables

Returns:
- plot bundle

## Full Pipeline

### `run_conseguiR()`

Purpose:
- run the full pipeline from input validation through final plot generation

Pipeline:
1. `validate_inputs()`
2. `prepare_germline_scores()`
3. `prepare_somatic_scores()`
4. `prepare_epigenomic_scores()`
5. `build_scored_gene_reg_graph()`
6. `run_gene_reg_diffusion()`
7. `plot_scores()` if the user wants intermediate score/diffusion plots
8. `call_selected_subgraph()`
9. `plot_selected_subgraph()`

Main outputs:
- all intermediate bundles
- final plot bundle

Returns:
- pipeline bundle

## Dependency Chain

The intended dependency chain is:

- `validate_inputs()` does not depend on downstream steps
- `prepare_germline_scores()`, `prepare_somatic_scores()`, and `prepare_epigenomic_scores()` depend only on validated raw inputs
- `build_scored_gene_reg_graph()` depends on the score bundles
- `run_gene_reg_diffusion()` depends on the scored graph bundle
- `plot_scores()` depends on one of the score/diffusion bundles or a tabular score object
- `call_selected_subgraph()` depends on the diffusion bundle
- `plot_selected_subgraph()` depends on the selected subgraph bundle
- `run_conseguiR()` orchestrates all of the above

## Design Notes

- Externals are task-oriented, not step-number-oriented
- Each layer accepts either explicit file paths or the previous layer's bundle
- Score preparation is intentionally split by method/modality rather than collapsed into one oversized function
- Downstream functions should not require raw upstream inputs once a bundle exists
