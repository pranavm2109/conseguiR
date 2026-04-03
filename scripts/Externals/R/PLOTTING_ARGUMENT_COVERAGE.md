# Plotting Argument Coverage

This document describes how `conseguiR` currently exposes arguments for the
selected-subgraph plotting and visualization-bundle stage.

The relevant user-facing functions are:

- `plot_selected_subgraph()`
- `run_conseguiR()`

Here, `bundle` refers to the returned visualization object with `objects`,
`output_paths`, and `config`, and `wrapper-level input` means a parameter
exposed directly by the user-facing R wrapper.

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
