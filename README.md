# conseguiR

`conseguiR` is an R package for integrative prioritization of genes and
regulatory elements using germline association signal, somatic mutation signal,
and epigenomic regulatory signal in a graph-based framework.

The package is designed for analyses where no single modality is sufficient on
its own. Instead of treating GWAS hits, somatic signals, and regulatory
activity as disconnected results, `conseguiR` maps them into a shared
gene-regulatory context, propagates signal through that context, and returns
interpretable locus-level and network-level outputs.

## What the package does

At a high level, `conseguiR` supports the following workflow:

1. score germline signal with MAGMA at genes and regulatory elements
2. score somatic signal with `dndscv` at genes and `fishHook` at regulatory
   elements
3. score epigenomic regulatory activity from bigWig tracks
4. impose all modality-specific scores onto a backend gene-regulatory graph
5. run diffusion to integrate regulatory evidence into gene-level support
6. call a compact selected subgraph on a gene-gene network
7. visualize scores, loci, and the final selected subgraph

The package returns structured R objects by default. Disk output is optional.

## Main exported workflow

The package exposes both stage-wise functions and a top-level wrapper.

### Core stage-wise functions

- `prepare_germline_scores()`
- `prepare_somatic_scores()`
- `prepare_epigenomic_scores()`
- `build_scored_gene_reg_graph()`
- `run_gene_reg_diffusion()`
- `call_selected_subgraph()`
- `plot_scores()`
- `plot_locus_context()`
- `plot_selected_subgraph()`

### End-to-end wrapper

- `run_conseguiR()`

The package also exports several lower-level branch-specific wrappers such as
`run_germline_gene_scoring()` and `run_somatic_regulatory_scoring()`. Those are
useful when a user wants direct branch-level control, but the functions above
are the main user-facing workflow.

## Quick start

### One-shot workflow

```r
result <- run_conseguiR(
  gwas_sumstats = "<path-or-table>",
  somatic_maf = "<path-or-table>",
  reference_bfile = "<plink-reference-prefix>",
  dndscv_refdb = "<dndscv-reference-db>",
  epigenomic_tracks = c("<track1.bw>", "<track2.bw>", "<track3.bw>"),
  target_genes = 50L,
  verbose = TRUE
)
```

This returns a pipeline bundle containing:

- germline score bundle
- somatic score bundle
- epigenomic score bundle
- scored gene-regulatory graph
- diffusion output
- selected subgraph
- final selected-subgraph plot bundle

By default, `conseguiR` uses its backend-owned ENCODE-derived regulatory
universe for regulatory scoring and graph imposition. In the common workflow,
you do not need to supply `reg_ref_path` manually.

### Stage-wise workflow

```r
germline <- prepare_germline_scores(
  gwas_sumstats = "<path-or-table>",
  reference_bfile = "<plink-reference-prefix>"
)

somatic <- prepare_somatic_scores(
  maf = "<path-or-table>",
  refdb = "<dndscv-reference-db>"
)

epigenomic <- prepare_epigenomic_scores(
  bw_files = c("<track1.bw>", "<track2.bw>", "<track3.bw>")
)

scored_graph <- build_scored_gene_reg_graph(
  germline_scores = germline,
  somatic_scores = somatic,
  epigenomic_scores = epigenomic
)

diffusion <- run_gene_reg_diffusion(scored_graph = scored_graph)

selected_subgraph <- call_selected_subgraph(
  diffusion = diffusion,
  target_genes = 50L
)

selected_plot <- plot_selected_subgraph(
  selected_subgraph = selected_subgraph,
  save_plot = FALSE
)
```

## Plotting

`conseguiR` supports three main plotting layers.

### 1. Score plots

Use `plot_scores()` for generic score plots.

```r
plot_scores(
  scores = germline,
  which = "gene_scores",
  plot_mode = "rank",
  save_plot = FALSE
)

plot_scores(
  scores = somatic,
  which = "gene_scores",
  plot_mode = "volcano",
  save_plot = FALSE
)
```

`plot_mode` controls geometry:

- `plot_mode = "rank"` for ranked feature plots
- `plot_mode = "volcano"` for z-score versus `-log10(p)` plots

### 2. Locus plots

Use `plot_locus_context()` to inspect a genomic region.

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
  save_plot = FALSE
)
```

This plot can show:

- GWAS locus SNPs
- regulatory-element input scores by modality
- regulatory-to-gene links in the locus
- post-diffusion gene support
- optional SNP labels and highlighted genes

### 3. Selected-subgraph plots

Use `plot_selected_subgraph()` for the final network view.

```r
plot_selected_subgraph(
  selected_subgraph = selected_subgraph,
  title = "Selected disease subgraph",
  save_plot = FALSE
)
```

All plotting functions return editable `ggplot` objects inside the returned
bundle.

## Package backend resources

`conseguiR` ships backend seed resources with the package under
`inst/extdata/backend`.

These resources support:

- the gene-regulatory graph
- the gene-gene graph
- gene and regulatory annotations needed by the workflow

In a source-checkout or development workflow, the package may materialize or
reuse working backend files in `data/processed`. In an installed-package
workflow, users interact with the package-owned backend resources rather than
rebuilding those graphs from scratch.

## External dependencies

`conseguiR` depends on a small number of nontrivial scientific tools.

### MAGMA

MAGMA is required for germline scoring and must be available externally.

The package looks for MAGMA in this order:

1. the explicit `magma_path` argument
2. `options(conseguiR.magma_path = "/path/to/magma")`
3. `Sys.getenv("CONSEGUIR_MAGMA_PATH")`
4. `magma` on `PATH`

These settings should point to the MAGMA executable itself, not just its
folder. The `PATH` fallback only works when MAGMA can already be invoked by
typing `magma` in the shell environment seen by R. This makes autodiscovery a
useful convenience on standard installs, but explicit configuration is often
safer on HPC or project-local setups.

### Python-backed stages

Diffusion and selected-subgraph calling use Python-backed stages managed
internally through `basilisk`. Users do not normally need to wire up a Python
interpreter by hand just to run the package, although developer and HPC
workflows may still use a project-specific conda environment for convenience.

## Documentation

For practical setup instructions, see:

- [INSTALL.md](INSTALL.md)

For the full package walkthrough, see:

- `vignettes/conseguiR-overview.Rmd`

The vignette demonstrates the intended stage-wise workflow and the top-level
`run_conseguiR()` workflow on real package-shaped inputs.

## Current package state

`conseguiR` is a research package under active development rather than a CRAN
release. The main public workflow is in place, the vignette is intended to be
reproducible, and the installed package is designed to ship the backend seeds
needed for normal use.

## Repository

- Source: <https://github.com/pranavm2109/conseguiR>
- Issues: <https://github.com/pranavm2109/conseguiR/issues>
