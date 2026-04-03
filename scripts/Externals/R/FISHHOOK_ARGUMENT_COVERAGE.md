# fishHook Argument Coverage

This document describes how `conseguiR` currently exposes `fishHook`
arguments through the external somatic regulatory-scoring API.

The relevant user-facing functions are:

- `run_somatic_regulatory_scoring()`
- `prepare_somatic_scores()`

Here, `wrapper-level input` means a parameter exposed directly by the
user-facing R wrapper, while `fishhook_args` is the passthrough list for
additional model options.

## Important Scope Note

`conseguiR` does **not** currently wrap every `fishHook` tuning parameter as a
named R argument.

Instead, `fishHook` arguments fall into three buckets:

1. explicitly supported as named wrapper arguments
2. supported through `fishhook_args`
3. not currently surfaced by the package workflow

This means users can still pass additional `fishHook` options through
`fishhook_args`, even when those options are not exposed as first-class named
arguments.

## Wrapper-Level Inputs

These arguments are exposed directly because they are central to the current
regulatory somatic workflow.

### `maf`

Somatic mutation input as a path or table.

This is harmonized internally into the mutation format used to build the
`GRanges` events supplied to `fishHook`.

### `reg_ref_path`

Regulatory-element reference path used to construct the hypothesis set for
`fishHook`.

### `eligible_gr`

Optional eligible territory as a `GRanges`.

If omitted, the current internal workflow constructs a default hg38 eligible
territory.

### `fishhook_covariates`

Optional list of `fishHook::Cov` objects or covariate specifications.

This is the direct path for advanced users who want to provide their own
`fishHook` covariates.

### `fishhook_covariate_data`

Optional tabular covariate data used by the helper path that derives default
`fishHook` covariates from regulatory-element metadata.

In the current workflow, this is the main route for supplying regulatory
covariates without hand-constructing `Cov` objects.

### `idcol`

Sample identifier column used in the `fishHook` call.

Current default:

- `Tumor_Sample_Barcode`

## Additional fishHook Arguments

### `fishhook_args`

Named list of additional arguments passed into the `fishHook` construction
path.

Use this when:

- you need a `fishHook` option that is not wrapped as a named external
  argument
- you want to tune the regulatory somatic model beyond the default wrapper
  surface

## How This Appears in the External API

### `run_somatic_regulatory_scoring()`

Exposes:

- `maf`
- `reg_ref_path`
- `eligible_gr`
- `fishhook_covariates`
- `fishhook_covariate_data`
- `idcol`
- `fishhook_args`

### `prepare_somatic_scores()`

Exposes the `fishHook` side of the joint somatic workflow through:

- `reg_ref_path`
- `eligible_gr`
- `fishhook_covariates`
- `fishhook_covariate_data`
- `fishhook_idcol`
- `fishhook_args`

These apply only to the regulatory `fishHook` branch of the somatic pipeline.

## What Is Not Yet Fully Wrapped

The package currently does **not** expose every possible `fishHook` option as a
named wrapper argument.

If a user needs an option not listed above, the intended path is:

1. pass it through `fishhook_args`
2. if it becomes commonly useful, promote it later to a named wrapper argument
