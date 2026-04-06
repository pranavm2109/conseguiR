# Plotting Argument Coverage

This document describes how `conseguiR` currently exposes arguments for:

- score plotting
- diffusion plotting
- selected-subgraph plotting and visualization bundles

The relevant user-facing functions are:

- `plot_scores()`
- `plot_diffusion()`
- `plot_selected_subgraph()`
- `run_conseguiR()`

Here, `bundle` refers to the returned visualization object with `objects`,
`output_paths`, and `config`, and `wrapper-level input` means a parameter
exposed directly by the user-facing R wrapper.

## Score Plotting

### `scores`

Optional score bundle returned by a scoring wrapper.

### `table`

Optional explicit in-memory score table.

### `which`

Optional score-table selector when the supplied bundle contains multiple score
tables.

Common values:

- `gene_scores`
- `reg_scores`

### `plot_type`

Supported values:

- `ranked_points`
- `top_bar`
- `histogram`

### `feature_column`

Optional explicit feature-label column. Leave `NULL` to let the plotting helper
infer the best available label column.

### `value_column`

Optional explicit numeric score column. Leave `NULL` to let the plotting helper
infer the score column.

### `highlight_features`

Optional character vector of features to highlight in the score plot.

## Diffusion Plotting

### `diffusion`

Optional diffusion bundle returned by `run_gene_reg_diffusion()`.

### `table`

Optional explicit in-memory diffusion table.

### `which`

Which diffusion table to plot.

Supported values:

- `all_genes`
- `top_genes`

### `gene_column`

Optional explicit gene-label column. Leave `NULL` to use the default/inferred
gene label column.

### `score_column`

Optional explicit numeric diffusion-score column.

### `highlight_genes`

Optional character vector of genes to highlight in the diffusion plot.

## Selected Subgraph Plotting

## Input Surface

### `selected_subgraph`

Optional selected-subgraph bundle returned by `call_selected_subgraph()`.

If provided, the wrapper resolves nodes, edges, and summary objects and paths
from the bundle automatically.

### `nodes`, `edges`, `summary`

Optional in-memory selected-subgraph tables.

### `nodes_path`, `edges_path`, `summary_path`

Optional explicit file paths for the selected-subgraph tables.

## Bundle-Saving Controls

### `bundle_output_prefix`

Output prefix used when saving the visualization bundle.

### `save_bundle`

Logical flag controlling whether the visualization bundle is saved.

## Figure-Saving Controls

### `plot_file_path`

Optional file path for the saved figure.

### `save_plot`

Logical flag controlling whether the figure is saved.

If `save_plot = TRUE`, `plot_file_path` must be provided.

## Plot Appearance Controls

### `title`

Custom plot title.

### `layout`

Graph layout name used for the selected-subgraph visualization.

### `top_n_labels`

Number of node labels to show.

Current default:

- `Inf`

## Output Rendering Controls

### `width`

Plot width in inches.

### `height`

Plot height in inches.

### `dpi`

Saved figure DPI.
