# conseguiR

`conseguiR` is an R package for integrating germline, somatic, and epigenomic
signals on gene-regulatory graphs, running diffusion, selecting a compact
gene-gene subgraph, and generating both summary and exploratory plots.

The package is designed to be object-first:

- compute functions return R objects and bundles by default
- file writing is optional and should happen only when the user explicitly asks
  for saved outputs
- plotting functions are the natural place to save figures

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
9. optionally inspect intermediate score outputs and specific genomic loci

The package ships prebuilt unscored backend graphs for:

- the gene-regulatory graph
- the gene-gene graph

These are treated as package-owned resources rather than user-provided inputs.

Diffusion and selected-subgraph calling are now run inside a managed Python
environment via `basilisk`. That means Python is still required for those two
stages, but the package no longer expects users to wire up a project-specific
interpreter by hand. The package provisions the required Python stack for those
steps internally.

Germline scoring still depends on MAGMA, but MAGMA is now treated as an
external prerequisite rather than a bundled package binary. You can supply it
explicitly with `magma_path`, register it once with
`options(conseguiR.magma_path = "/path/to/magma")`, set
`CONSEGUIR_MAGMA_PATH`, or make `magma` available on your system `PATH`.

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
- `call_selected_subgraph()`
- `plot_locus_context()`
- `plot_selected_subgraph()`
- `run_conseguiR()`

Detailed API notes are in:

- [scripts/Externals/R/EXTERNAL_API_OVERVIEW.md](scripts/Externals/R/EXTERNAL_API_OVERVIEW.md)

## Quickstart

The highest-level workflow is through `run_conseguiR()`:

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

By default, this returns a master pipeline bundle in memory. If you want the
pipeline to also save artifacts to disk, pass an explicit `output_dir`.

If you want a staged workflow instead of the one-shot wrapper, the intended
sequence is:

```r
germline <- prepare_germline_scores(...)
somatic <- prepare_somatic_scores(...)
epigenomic <- prepare_epigenomic_scores(...)

scored_graph <- build_scored_gene_reg_graph(
  germline_scores = germline,
  somatic_scores = somatic,
  epigenomic_scores = epigenomic
)

diffusion <- run_gene_reg_diffusion(scored_graph = scored_graph)
selected_subgraph <- call_selected_subgraph(diffusion = diffusion)
```

## Plotting

`conseguiR` currently supports three complementary plotting modes:

1. score-stage plots
2. selected-subgraph plots
3. exploratory locus plots

Intermediate score plots:

```r
plot_germline_gene_scores(germline_scores = germline, stage = "pre")
plot_germline_gene_scores(
  germline_scores = germline,
  diffusion = diffusion,
  stage = "post"
)

plot_somatic_gene_scores(somatic_scores = somatic, stage = "pre")
plot_somatic_gene_scores(
  somatic_scores = somatic,
  diffusion = diffusion,
  stage = "post"
)

plot_epigenomic_gene_scores(diffusion = diffusion)
plot_germline_reg_scores(germline_scores = germline)
plot_somatic_reg_scores(somatic_scores = somatic)
plot_epigenomic_reg_scores(epigenomic_scores = epigenomic)
```

Selected-subgraph plot:

```r
plot_selected_subgraph(
  selected_subgraph = selected_subgraph,
  plot_file_path = "selected_subgraph.pdf"
)
```

Exploratory locus plot:

```r
plot_locus_context(
  chromosome = "8",
  start = 127200000,
  end = 128200000,
  scored_graph = scored_graph,
  diffusion = diffusion,
  selected_subgraph = selected_subgraph,
  gwas_sumstats = "<path-or-table>",
  label_top_lit_snps = 3L,
  plot_file_path = "MYC_locus_context.pdf"
)
```

`plot_locus_context()` shows:

- regulatory-element somatic, epigenomic, and germline z-score tracks
- a combined regulatory score track
- a post-diffusion gene score track
- regulatory-element to gene links within the locus
- optional SNP labels, preferring literature-backed SNPs when available and
  falling back to top GWAS SNPs in the top germline regulatory elements

## Installation

For installation and setup notes, see:

- [INSTALL](INSTALL)

## Current status

This is a strong research-grade draft rather than a polished release.

The main remaining rough edges are:

- some package-facing wrappers still delegate into `scripts/...`
- MAGMA still matters for germline scoring
- the `basilisk`-managed Python path for
  `run_gene_reg_diffusion()` and `call_selected_subgraph()` has been validated
  in this development environment and still deserves one clean-machine check
  before release
- the package still needs some final polish around fresh-clone validation and
  resource discovery
