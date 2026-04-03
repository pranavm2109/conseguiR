# Subgraph Argument Coverage

This document describes how `conseguiR` currently exposes arguments for the
Python-backed subgraph-calling stage.

The relevant user-facing functions are:

- `call_selected_subgraph()`
- `run_conseguiR()`

Here, `wrapper-level input` means a parameter exposed directly by the
user-facing R wrapper rather than hidden inside the internal Python module.

## Important Scope Note

The subgraph-calling stage is also exposed almost entirely through named wrapper
arguments.

So this document describes the solver and candidate-selection surface that the
wrapper already makes available directly.

## Input Surface

### `diffusion`

Optional diffusion bundle returned by `run_gene_reg_diffusion()`.

If provided, the wrapper resolves the diffusion table path automatically.

### `diffusion_path`

Explicit diffusion-results path.

### `gg_nodes_path`

Gene-gene node table path.

### `gg_edges_path`

Gene-gene edge table path.

## Output Controls

### `output_dir`

Directory where selected-subgraph outputs will be written.

### `output_stem`

Stem used to construct the selected-subgraph output files.

## Cardinality and Candidate Controls

### `target_genes`

Requested subgraph cardinality.

### `candidate_pool_size`

Number of top candidate genes passed into the optimization problem before edge
restriction and exact cardinality solving.

### `min_confidence`

Minimum gene-gene edge confidence allowed into the solver model.

### `max_edges_in_model`

Upper bound on the number of candidate edges retained in the optimization
problem.

## Objective Controls

### `node_prize_weight`

Weight applied to node prizes in the optimization objective.

### `edge_conf_weight`

Weight applied to edge-confidence reward in the optimization objective.

### `edge_cost_weight`

Weight applied to edge-cost penalty in the optimization objective.

### `node_scale`

Integer scaling factor applied to node prizes before CP-SAT optimization.

### `edge_scale`

Integer scaling factor applied to edge objective terms before CP-SAT
optimization.

## Solver Controls

### `max_time_seconds`

Time limit for the CP-SAT solve.

### `num_workers`

Number of solver workers.

### `random_seed`

Solver random seed.

## Column-Mapping Controls

### `prize_column`

Diffusion column used as the node prize.

### `confidence_column`

Gene-gene edge column interpreted as connection confidence.

### `edge_cost_column`

Gene-gene edge column interpreted as edge cost or penalty.

## Runtime Control

### `python_path`

Optional explicit Python interpreter path.
