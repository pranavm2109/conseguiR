# Installation

This document describes the intended installation path for `conseguiR` and the
small number of external requirements that matter for a successful first run.

## What you need

`conseguiR` combines:

- R
- package-managed Python-backed stages via `basilisk`
- MAGMA for germline scoring
- package-shipped backend graph resources

For most users, the two main setup tasks are:

1. install the R package and its dependencies
2. make MAGMA available

## 1. Install the package

### From GitHub

```r
remotes::install_github("pranavm2109/conseguiR", dependencies = TRUE)
```

Then load it with:

```r
library(conseguiR)
```

### From a local clone

```r
devtools::install()
```

or for active development:

```r
devtools::load_all()
```

## 2. MAGMA

MAGMA is required for germline scoring.

The package looks for MAGMA in the following order:

1. the explicit `magma_path` argument passed to a germline scoring function
2. `options(conseguiR.magma_path = "/path/to/magma")`
3. `Sys.getenv("CONSEGUIR_MAGMA_PATH")`
4. `magma` on your system `PATH`

Important:

- these settings should point to the MAGMA executable itself, not just the
  containing folder
- for example, use `/path/to/magma` rather than `/path/to/`
- the `PATH` fallback only works if MAGMA can be run directly by typing
  `magma` in the same shell environment seen by R

In other words, autodiscovery is mainly a convenience for systems where MAGMA
has already been installed in a standard executable location or has already
been added to `PATH`. On HPC systems or project-specific installs, it is often
safer to set the executable path explicitly.

Examples:

```r
options(conseguiR.magma_path = "/path/to/magma")
```

or in the shell:

```bash
export CONSEGUIR_MAGMA_PATH=/path/to/magma
```

For example, on HPC:

```bash
export CONSEGUIR_MAGMA_PATH="/path/to/project/tools/magma_v1/magma"
```

If MAGMA is unavailable, germline scoring functions should fail with a clear
message rather than silently misbehaving.

## 3. Python-backed stages

`run_gene_reg_diffusion()` and `call_selected_subgraph()` use Python-backed
stages managed internally through `basilisk`.

In a normal installed-package workflow, users do not usually need to configure
Python manually. The package handles the managed Python environment for those
stages.

For developer or HPC workflows, a project-specific conda environment can still
be helpful for reproducibility, but it is not the main package interface.

## 4. Backend graph resources

`conseguiR` ships backend seed resources with the package under:

- `inst/extdata/backend`

These packaged resources support the installed-package workflow and include the
backend graph seeds needed for normal use.

In a source-checkout or development session, the package may also materialize
working backend files in:

- `data/processed`

That working directory is for local or development use. Installed users should
think of `inst/extdata/backend` as the package-owned backend resource layer.

## 5. First-run sanity check

After installation, a clean first check is:

```r
library(conseguiR)
check_conseguiR_runtime()
initialize_backend_graphs()
```

What to expect:

- `library(conseguiR)` should load cleanly
- `check_conseguiR_runtime()` should report a healthy runtime
- `initialize_backend_graphs()` should seed or reuse backend resources

If you want to inspect shipped backend resources directly:

```r
list.files(system.file("extdata", "backend", package = "conseguiR"))
```

## 6. Minimal usage check

Once the package loads successfully, a minimal end-to-end call looks like:

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

This returns the stage bundles, selected subgraph, and final graph plot bundle
in memory. If you also want saved outputs, provide `output_dir`.

By default, `conseguiR` uses its backend-owned ENCODE-derived regulatory
universe for regulatory scoring and graph imposition. In the common workflow,
you do not need to supply `reg_ref_path` manually.

## 7. Vignette and documentation

The current public documentation is centered on:

- this installation guide
- function-level help pages
- the package README

A refreshed installed-package vignette will be added back after the release
cleanup pass.

## 8. Troubleshooting notes

### If germline scoring fails

Check:

- that MAGMA is installed and executable
- that `reference_bfile` points to a valid PLINK prefix
- that the GWAS summary-statistics inputs match what MAGMA expects

### If diffusion or subgraph calling fails

Check:

- that `check_conseguiR_runtime()` reports a healthy runtime
- that the package was loaded in a clean session
- that backend graph resources were initialized successfully

### If you are working from a source checkout

Remember that source-mode sessions may use a local working backend in
`data/processed`, while installed-package sessions rely on the shipped package
resources under `inst/extdata/backend`.

## 9. Recommended installation story

The intended user story is:

1. install the package from GitHub or from a local source bundle
2. make MAGMA available
3. load the package
4. run `check_conseguiR_runtime()`
5. call `run_conseguiR()` directly or work stage-wise through the exported API

That is the level of setup the package is aiming to support reliably.
