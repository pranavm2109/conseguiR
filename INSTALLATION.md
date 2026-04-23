# Installation

This document describes the current installation story for `conseguiR`.

It is intentionally practical and honest: the package is usable, but the setup
is not yet fully polished into a one-command install.

## Current status

At the moment, `conseguiR` is best thought of as:

- an R package with exported user-facing wrappers
- a mixed R + Python workflow
- a package that also depends on external scientific tools and resources

So installation currently involves four layers:

1. R
2. Python
3. MAGMA
4. package-owned backend graph resources

## 1. Clone the repository

```bash
git clone <your-repo-url>
cd conseguiR
```

## 2. Install the R package dependencies

At minimum, the package metadata currently declares these core R dependencies:

- `data.table`
- `dplyr`
- `jsonlite`
- `igraph`
- `ggplot2`
- `GenomicRanges`
- `IRanges`
- `GenomeInfoDb`
- `AnnotationDbi`
- `S4Vectors`

And the workflow also relies on additional packages used in the scoring and
plotting stages, including:

- `testthat`
- `reticulate`
- `ggrepel`
- `tidygraph` (optional)
- `RCy3` (optional, for Cytoscape export)
- `rtracklayer`
- `BSgenome.Hsapiens.UCSC.hg38`
- `dndscv`
- `fishHook`

Depending on how you prefer to work, you can install the package in
development mode with `devtools`:

```r
devtools::load_all()
```

or install it locally:

```r
devtools::install()
```

## 3. Prepare the Python environment

The recommended Python environment name is:

- `lymphoma_graph_env`

The repository now includes a matching conda specification:

- `environment.yml`

So on a fresh machine or HPC node, the intended first step is:

```bash
conda env create -f environment.yml
conda activate lymphoma_graph_env
```

This environment is used primarily for:

- diffusion
- subgraph calling
- developer and HPC execution of helper scripts such as backend-seed builders

So in the current workflow, the easiest path is to ensure that
`lymphoma_graph_env` exists, is active for your shell session, and is visible
to `conda`.

## 4. Ensure MAGMA is available

The germline scoring stage uses MAGMA.

`conseguiR` resolves MAGMA in this order:

- the explicit `magma_path` argument passed to a germline scoring function
- `options(conseguiR.magma_path = "/path/to/magma")`
- `Sys.getenv("CONSEGUIR_MAGMA_PATH")`
- `magma` on your system `PATH`

So before running germline scoring, make sure that:

- a MAGMA 1.1 binary is available through one of those paths
- it is executable
- the related reference inputs you plan to use are available

## 5. Understand the current backend-resource expectation

`conseguiR` works with backend graph resources, especially:

- a gene-regulatory graph without user-specific scores
- a gene-gene graph used after diffusion

These backend graph resources are now intended to ship with the package as
package-owned resources.

On package load, `conseguiR` attempts to ensure these backend graph files exist
in the working backend directory. If they are missing there, it first tries to
seed them from the packaged backend graph resources.

At the moment, the gene-gene packaged seed is ready and the ENCODE gene-reg
compact seed is generated as a one-time developer step before distribution.

## 6. Sanity check package startup

After loading the package, a useful first check is:

```r
library(conseguiR)
check_conseguiR_runtime()
```

If the Python setup is healthy, the returned runtime status should report the
managed Python stage as available.

## 7. Minimal workflow check

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

The final plotted graph path is then available from:

- `result$plot$output_paths$plot_file_path`

## 8. Fresh-clone validation checklist

After cloning the repository into a separate directory, a good first validation
path is:

1. install the package with `devtools::install()`
2. load the package with `library(conseguiR)`
3. run `check_conseguiR_runtime()`
4. run `initialize_backend_graphs()`

In a clean R session, that looks like:

```r
devtools::install()
library(conseguiR)

check_conseguiR_runtime()
initialize_backend_graphs()
```

The ideal success state is:

- `library(conseguiR)` loads cleanly
- `check_conseguiR_runtime()` reports `OK`
- `initialize_backend_graphs()` reports `reused` or `seeded`

Useful follow-up checks:

```r
check_conseguiR_runtime(quiet = TRUE)
list.files(system.file("extdata", "backend", package = "conseguiR"))
```

If you want to go one step further, then try one exported workflow function
with real study inputs, for example `run_conseguiR(...)`.

## 9. What is still rough right now

The main installation-related caveats right now are:

- some package-facing wrappers still delegate into `scripts/...`
- Python setup is important for a smooth first run
- MAGMA is an external dependency rather than a pure R package dependency
- installation and package-resource discovery still need final cleanup

## Recommended next improvement

The next installation improvement for the package should be:

- a cleaner fully packaged install path with less reliance on development-time
  repository structure

That would move the package closer to the ideal user experience where loading
the library and calling the exported wrappers is enough, without relying on any
personal development-time folder layout.
