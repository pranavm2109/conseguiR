# dndscv Argument Coverage

This document describes how `conseguiR` currently exposes `dndscv` arguments
through the external somatic gene-scoring API.

The relevant user-facing functions are:

- `run_somatic_gene_scoring()`
- `prepare_somatic_scores()`

Here, `wrapper-level input` means a parameter exposed directly by the
user-facing R wrapper, while `dndscv_args` is the passthrough list for
additional model options.

## Important Scope Note

`conseguiR` does **not** currently wrap every `dndscv` parameter as a named R
argument.

Instead, `dndscv` arguments fall into three buckets:

1. explicitly supported as named wrapper arguments
2. supported through `dndscv_args`
3. not currently surfaced by the package workflow

This means users can still pass additional `dndscv` options through
`dndscv_args`, even when those options are not exposed as first-class named
arguments.

## Wrapper-Level Inputs

These arguments are exposed directly by the package wrappers because they are
central to the current somatic gene-scoring workflow.

### `maf`

Somatic mutation input as a path or table.

This is harmonized internally into the `dndscv` mutation format before model
fitting.

### `refdb`

Reference database path passed to `dndscv`.

This is required in the current workflow.

The package also uses `refdb` to detect chromosome naming style and harmonize
the input MAF accordingly before the `dndscv` call.

### `cv`

Optional `dndscv` covariate input.

This is passed through directly as the `cv` argument to `dndscv`.

### `max_muts_per_gene_per_sample`

Wrapper-level control for the `dndscv` per-gene mutation cap.

Passed directly to:

- `max_muts_per_gene_per_sample`

### `max_coding_muts_per_sample`

Wrapper-level control for the `dndscv` per-sample coding mutation cap.

Passed directly to:

- `max_coding_muts_per_sample`

## Additional dndscv Arguments

### `dndscv_args`

Named list of additional arguments passed through to `dndscv` via `do.call()`.

Use this when:

- you need a `dndscv` parameter that is not wrapped as a named external
  argument
- you want to override or extend the default gene-scoring call

Important default:

- unless the user overrides it explicitly, `conseguiR` requests
  `onesided = TRUE` when the installed `dndscv` version supports that
  argument
- when one-sided output columns such as `ppos_cv` and `pneg_cv` are returned,
  `conseguiR` uses those directional p-values directly for z-score extraction

## How This Appears in the External API

### `run_somatic_gene_scoring()`

Exposes:

- `maf`
- `refdb`
- `cv`
- `max_muts_per_gene_per_sample`
- `max_coding_muts_per_sample`
- `dndscv_args`

### `prepare_somatic_scores()`

Exposes the `dndscv` side of the joint somatic workflow through:

- `gene_cv`
- `gene_max_muts_per_gene_per_sample`
- `gene_max_coding_muts_per_sample`
- `dndscv_args`

These apply only to the gene-level `dndscv` branch of the somatic pipeline.

## What Is Not Yet Fully Wrapped

The package currently does **not** expose every possible `dndscv` option as a
named wrapper argument.

If a user needs an option not listed above, the intended path is:

1. pass it through `dndscv_args`
2. if it becomes commonly useful, promote it later to a named wrapper argument
