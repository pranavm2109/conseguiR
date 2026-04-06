# conseguiR

`conseguiR` is an R package for integrating germline, somatic, and epigenomic
signals on gene-regulatory graphs, running diffusion, selecting a compact
gene-gene subgraph, and exporting a final visualization bundle.

This repository contains a first-draft end-to-end implementation of that
workflow, including exported user-facing wrappers, internal pipeline stages,
and plotting support.

## Overview

The current package workflow is:

1. validate user inputs
2. compute germline scores with MAGMA
3. compute somatic scores with `dndscv` and `fishHook`
4. compute regulatory epigenomic scores from bigWig tracks
5. impose all scores onto the backend gene-regulatory graph
6. run diffusion on the scored gene-regulatory graph
7. call a cardinality-constrained gene-gene subgraph
8. build a visualization bundle and save a plotted graph

The package ships prebuilt unscored backend graphs for:

- the gene-regulatory graph
- the gene-gene graph

These are treated as package-owned resources rather than user-provided inputs.

## User API

The package-facing API currently exports:

- `validate_inputs()`
- `run_germline_gene_scoring()`
- `run_germline_regulatory_scoring()`
- `prepare_germline_scores()`
- `run_somatic_gene_scoring()`
- `run_somatic_regulatory_scoring()`
- `prepare_somatic_scores()`
- `prepare_epigenomic_scores()`
- `build_scored_gene_reg_graph()`
- `run_gene_reg_diffusion()`
- `plot_scores()`
- `plot_diffusion()`
- `call_selected_subgraph()`
- `plot_selected_subgraph()`
- `run_conseguiR()`

Detailed API notes are in:

- [scripts/Externals/R/EXTERNAL_API_OVERVIEW.md](scripts/Externals/R/EXTERNAL_API_OVERVIEW.md)

## Quickstart

The intended public workflow is through the exported wrappers, for example:

```r
result <- run_conseguiR(
  gwas_sumstats = "<path-or-table>",
  somatic_maf = "<path-or-table>",
  reg_ref_path = "<regulatory-reference-path>",
  reference_bfile = "<plink-reference-prefix>",
  dndscv_refdb = "<dndscv-refdb-path>",
  epigenomic_track_dir = "<bigwig-directory>",
  target_genes = 50L
)
```

The final plot path is then available from:

- `result$plot$output_paths$plot_file_path`

## Installation

For installation and setup notes, see:

- [INSTALLATION.md](INSTALLATION.md)

## Current status

This is a strong first draft rather than a polished release.

The main remaining rough edges are:

- some package-facing wrappers still delegate into `scripts/...`
- Python and MAGMA setup still matter for a smooth first run
- the package still needs some final packaging polish around installation and
  resource discovery
