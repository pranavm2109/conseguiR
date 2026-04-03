# Diffusion Argument Coverage

This document describes how `conseguiR` currently exposes arguments for the
Python-backed diffusion stage.

The relevant user-facing functions are:

- `run_gene_reg_diffusion()`
- `run_conseguiR()`

Here, `wrapper-level input` means a parameter exposed directly by the
user-facing R wrapper rather than hidden inside the internal Python module.

## Important Scope Note

Unlike MAGMA, `dndscv`, or `fishHook`, the diffusion stage is already exposed
almost entirely as named wrapper arguments.

So this document is mainly a map of what the wrapper controls rather than a
split between named arguments and an `extra_args` passthrough.

## Input Surface

### `scored_graph`

Optional scored graph bundle returned by `build_scored_gene_reg_graph()`.

If provided, the wrapper resolves the scored node and edge paths from the
bundle automatically.

### `nodes_path`

Explicit scored node table path.

### `edges_path`

Explicit scored edge table path.

Use `nodes_path` and `edges_path` when you want to run diffusion without
passing the full scored-graph bundle.

## Output Controls

### `output_dir`

Directory where diffusion outputs will be written.

### `output_stem`

Stem used to construct:

- `<output_stem>_all_genes.tsv`
- `<output_stem>_topN.tsv`

## Diffusion Hyperparameters

### `top_k`

Number of top-confidence regulatory neighbors considered per gene during the
controlled regulatory-to-gene smoothing step.

### `confidence_power`

Exponent applied to edge confidence before weighting incoming regulatory signal.

### `beta_germline`

Weight applied to incoming germline regulatory signal before combining it with
the gene-level germline score.

### `beta_somatic`

Weight applied to incoming somatic regulatory signal before combining it with
the gene-level somatic score.

### `beta_epigenomic`

Weight applied to incoming epigenomic regulatory signal before combining it
with the gene-level epigenomic score.

### `positive_only`

Logical flag controlling whether only positive incoming regulatory signal is
allowed to contribute to the diffusion aggregate.

### `reg_signal_clip`

Maximum absolute regulatory signal used before aggregation.

This limits the influence of extreme regulatory scores.

### `top_n_to_save`

Number of top-ranked genes written to the separate top-gene diffusion output.

## Runtime Control

### `python_path`

Optional explicit Python interpreter path.

If omitted, the wrapper will use the package’s Python-resolution logic.
