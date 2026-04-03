# Installation

This document describes the current first-draft installation story for
`conseguiR`.

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
4. backend graph/resources expected by the package workflow

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

This environment is used primarily for:

- diffusion
- subgraph calling

The package startup logic tries to discover the Python interpreter for this
environment and record it in:

- `getOption("conseguiR.python")`

So in the current first draft, the easiest path is to ensure that
`lymphoma_graph_env` exists and is visible to `conda`.

## 4. Ensure MAGMA is available

The germline scoring stage uses MAGMA.

The current internal defaults expect the executable at:

- `tools/magma_v1/magma`

So before running germline scoring, make sure that:

- the MAGMA binary is present
- it is executable
- the related reference inputs you plan to use are available

## 5. Understand the current backend-resource expectation

`conseguiR` works with backend graph resources, especially:

- a gene-regulatory graph without user-specific scores
- a gene-gene graph used after diffusion

The intended package behavior is that these backend graph resources should be
package-managed rather than tied to a private development workflow.

For this first draft, that story is still in progress. So if you are trying the
package before that final cleanup is complete, make sure the required backend
graph resources exist where the wrappers expect them.

## 6. Sanity check package startup

After loading the package, a useful first check is:

```r
library(conseguiR)
getOption("conseguiR.python")
```

If the Python setup is healthy, that option should resolve to a real Python
binary.

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

## 8. What is still rough in the first draft

The main installation-related caveats right now are:

- some package-facing wrappers still delegate into `scripts/...`
- Python setup is important for a smooth first run
- MAGMA is an external dependency rather than a pure R package dependency
- backend graph/resource management is not yet as automated as it should be

## Recommended next improvement

The next installation improvement for the package should be:

- a cleaner setup path for backend graph/resource creation and discovery

That would move the package closer to the ideal user experience where loading
the library and calling the exported wrappers is enough, without relying on any
personal development-time folder layout.
