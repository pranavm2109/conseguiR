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

What it is and is not:
- this is a lightweight sanity check
- it is meant to catch bad file paths, unreadable files, and obviously wrong
  columns before scoring starts
- it is not meant to fully dry-run MAGMA, dndscv, fishHook, or diffusion

Practical minimums:
- GWAS: one SNP identifier column, one chromosome column, one base-pair
  position column, and one p-value column
- somatic MAF: sample ID, chromosome, start, end, ref, alt
- regulatory reference: at least regulatory-element ID, chromosome, start, end
- epigenomic: at least three readable bigWig tracks

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

How to think about the two MAGMA stages:
- step 1 decides which SNPs belong to each feature
- step 2 decides how the annotated SNP-level signal is collapsed into one
  feature-level statistic

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

Common interpretation:
- use `sample_size` when the full GWAS has one fixed N
- use `sample_size_col` when N varies per row and the GWAS already contains
  that information
- in a typical run you use one or the other rather than supplying both
- `pval_duplicate` controls how MAGMA handles duplicate SNP IDs
- `reference_bfile` should be the shared PLINK prefix without file suffixes,
  e.g. `/path/to/g1000_eur/g1000_eur`
- `gene_model = "snp-wise=mean"` and `pval_use = c("SNP", "P")` are practical
  starting settings for many GWAS tables

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

Practical difference from the gene run:
- the "features" here are regulatory elements rather than genes
- users often keep the regulatory `annotation_window` narrower than the gene
  one because REs are already localized intervals

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

Practical decision rules:
- use the same LD reference (`reference_bfile`) for both runs unless you have a
  very specific reason not to
- keep the regulatory `annotation_window` narrower when regulatory elements are
  already localized intervals
- use `shared_args` only for wrapper-level settings you genuinely want to send
  to both MAGMA branches

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

Useful interpretation:
- this is the gene-centric somatic branch
- dndscv is where coding mutation burden and selection-like signal are
  summarized at gene level
- the most common extra argument users tune is `sm`
- if `cv` is supplied, it should already be in a dndscv-ready format; the
  package does not reshape arbitrary covariate tables for dndscv automatically
- `refdb` must be a dndscv-compatible reference `.rda` built for the same
  genome build as the MAF being analyzed

Minimal examples:
- `dndscv_args = list(sm = "192r_3w", kc = "cgc81")`
- `cv = NULL`

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

Useful interpretation:
- this is the regulatory-element-centric somatic branch
- `eligible_gr` defines the territory fishHook treats as eligible background
- `fishhook_covariate_data` is where users supply one row per regulatory
  element with covariates such as accessibility, replication, or GC-related
  quantities if they have them

Practical formatting:
- `fishhook_covariate_data` should be one row per regulatory element
- it should include a regulatory-element identifier column plus the covariates
  you want fishHook to model
- `fishhook_covariates` should already be a fishHook-ready specification if
  you provide it
- the regulatory-element IDs in `fishhook_covariate_data` should correspond to
  the IDs in `reg_ref_path`
- `idcol` should match the sample identifier field used by the somatic table
  after harmonization, usually `Tumor_Sample_Barcode`

Minimal example:
- `covariate_dt = data.frame(reg_elem_id = c("EH38E0080197", "EH38E2084302"), accessibility = c(1.2, 0.4), replication_timing = c(0.7, -0.1), gc_content = c(0.44, 0.51))`

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

Mental model:
- `prepare_somatic_scores()` does not fit one combined somatic model
- it runs dndscv and fishHook separately, then returns both outputs together

Formatting note:
- the same covariate rules apply here as in the lower-level dndscv and
  fishHook wrappers

## Epigenomic Scoring

### `prepare_epigenomic_scores()`

Purpose:
- compute regulatory-element epigenomic z scores from bigWig signal tracks

Main inputs:
- `reg_ref_path`
- `bw_files` or `track_dir`
- `min_tracks`
- `transform`

Main outputs:
- regulatory epigenomic score table
- optional diagnostics

Returns:
- epigenomic score bundle

Important semantic point:
- the current epigenomic score is a cross-track variability score, not a
  generic activity score
- high score means the regulatory element varies strongly across the supplied
  tracks and is therefore treated as more context-specific

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

Current default behavior:
- missing modalities are zero-filled rather than causing the feature to be
  dropped from the graph before diffusion
- this keeps the scored graph broader and avoids discarding biologically real
  genes or regulatory elements just because one modality is missing

Explicit table formats:
- gene score tables are safest when they look like `gene_id`, `zstat`
- regulatory score tables are safest when they look like `reg_elem_id`, `zstat`

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

How to think about the main knobs:
- `beta_germline`, `beta_somatic`, and `beta_epigenomic` control modality
  weighting
- `confidence_power` controls how strongly edge confidence affects propagation
- `top_k` controls how many top regulatory contributors influence each gene
- `reg_signal_clip` caps extreme regulatory signal before propagation

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

How to think about the main knobs:
- `target_genes` is the size of the final selected subgraph
- `candidate_pool_size` is the larger pool considered before the final
  selection
- `node_prize_weight` rewards strong diffusion signal
- `edge_conf_weight` rewards confident gene-gene edges
- `edge_cost_weight` penalizes expensive edges

Column-format examples:
- `prize_column = "post_integrated"`
- `confidence_column = "confidence"`
- `edge_cost_column = "weight"`

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

### `plot_locus_context()`

Purpose:
- create an exploratory locus-centered panel that combines regulatory-element
  z-score tracks, a combined regulatory score track, post-diffusion gene
  scores, regulatory links, and optional SNP labels

Main inputs:
- a merged post-diffusion gene-reg graph object, or a scored graph plus
  diffusion bundle/path
- locus coordinates (`chromosome`, `start`, `end`)
- optional selected-subgraph inputs for highlighting/filtering
- optional GWAS summary statistics for SNP labeling
- optional `rsid_pmid` table or a disease term such as `"DLBCL"` or
  `"lymphoma"` for LitVar-backed SNP labeling

Behavior:
- the top three tracks show regulatory-element somatic, epigenomic, and
  germline z-scores
- the `Reg elements` track colors regulatory elements by their combined
  pre-diffusion norm
- the bottom gene track colors genes by post-diffusion `conseguiR` score
- SNP labels prefer literature-backed SNPs and otherwise fall back to top GWAS
  SNPs in the top germline regulatory elements

Main outputs:
- ggplot object
- locus plotting bundle with feature/link tables
- optional SNP-label metadata and optional saved plot file

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
10. `plot_locus_context()` if the user wants an exploratory locus-level panel

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
- `plot_locus_context()` depends on a post-diffusion gene-reg graph view, plus
  optional GWAS/literature inputs for SNP labeling
- `run_conseguiR()` orchestrates all of the above

## Design Notes

- Externals are task-oriented, not step-number-oriented
- Each layer accepts either explicit file paths or the previous layer's bundle
- Score preparation is intentionally split by method/modality rather than collapsed into one oversized function
- Downstream functions should not require raw upstream inputs once a bundle exists
