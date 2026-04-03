# MAGMA Argument Coverage

This document describes how `conseguiR` currently exposes MAGMA arguments
through the external germline-scoring API.

Here, `wrapper-level input` means a parameter exposed directly by the
user-facing R wrapper, while `extra_args` is the passthrough vector for
additional MAGMA CLI options.

It is intentionally organized around the two MAGMA stages used by the package:

1. step 1: annotation
2. step 2: gene analysis

For each stage, there are two separate user-facing parameter bundles:

- gene run
- regulatory-element run

In `prepare_germline_scores()`, this becomes four argument groups:

- `gene_step1_args`
- `gene_step2_args`
- `reg_step1_args`
- `reg_step2_args`

## Important Scope Note

`conseguiR` does **not** currently wrap every MAGMA manual option as a named R
argument.

Instead, MAGMA arguments fall into three buckets:

1. explicitly supported as named R arguments
2. supported through `extra_args`
3. not currently surfaced by the package workflow

This means users can still reach additional MAGMA options through `extra_args`,
even when those options are not exposed as first-class named R parameters.

## Step 1: Annotation

Step 1 is driven by `run_magma_step1_annotation()` internally and exposed
through:

- `run_germline_gene_scoring(..., step1_args = list(...))`
- `run_germline_regulatory_scoring(..., step1_args = list(...))`
- `prepare_germline_scores(..., gene_step1_args = list(...), reg_step1_args = list(...))`

### Named Step 1 Arguments

#### `annotation_window`

Controls MAGMA annotation window behavior.

Expected forms:

- `NULL`
- a scalar window size
- a length-2 vector for asymmetric upstream/downstream windows

Mapped to MAGMA:

- `--annotate window=...`

#### `filter_path`

Path to a MAGMA-compatible filter file used during annotation.

Mapped to MAGMA:

- `--annotate filter=...`

#### `ignore_strand`

Logical flag controlling whether strand should be ignored during annotation.

Mapped to MAGMA:

- `--annotate ignore-strand`

#### `nonhuman`

Logical flag for non-human annotation mode.

Mapped to MAGMA:

- `--annotate nonhuman`

#### `extra_args`

Character vector of additional raw MAGMA step 1 arguments appended to the
annotation command.

Use this when:

- you need a MAGMA step 1 option that is not wrapped as a named R argument
- you want to pass through an uncommon manual option directly

## Step 2: Gene Analysis

Step 2 is driven by `run_magma_step2_gene_analysis()` internally and exposed
through:

- `run_germline_gene_scoring(..., step2_args = list(...))`
- `run_germline_regulatory_scoring(..., step2_args = list(...))`
- `prepare_germline_scores(..., gene_step2_args = list(...), reg_step2_args = list(...))`

### Core Step 2 Inputs Exposed at the Wrapper Level

These are not inside `step2_args`; they are direct wrapper arguments because
they are central to the pipeline setup.

#### `sample_size`

Fixed sample size for MAGMA gene analysis.

Mapped to MAGMA `--pval` syntax via:

- `N=<value>`

#### `sample_size_col`

Column name containing per-row sample sizes.

Mapped to MAGMA `--pval` syntax via:

- `N=<column_name>`

Only one of `sample_size` and `sample_size_col` should be supplied.

### Named Step 2 Arguments

#### `gene_model`

Model specification for MAGMA gene analysis.

Mapped to MAGMA:

- `--gene-model <value>`

#### `genes_only`

Logical flag controlling whether `--genes-only` is added.

This is commonly used in the current package workflow.

Mapped to MAGMA:

- `--genes-only`

#### `pval_use`

Character vector describing which columns MAGMA should use from the p-value
file.

Typical default:

- `c("SNP", "P")`

Mapped to MAGMA:

- `--pval ... use=SNP,P`

#### `pval_duplicate`

Optional duplicate-handling mode for the p-value file.

Mapped to MAGMA:

- `--pval ... duplicate=<value>`

#### `bfile_synonyms`

Optional path to a SNP synonym file for the reference `bfile`.

Mapped to MAGMA:

- `--bfile ... synonyms=<path>`

#### `bfile_synonym_dup`

Optional duplicate-handling mode for synonym matching.

Mapped to MAGMA:

- `--bfile ... synonym-dup=<value>`

#### `extra_args`

Character vector of additional raw MAGMA step 2 arguments appended to the gene
analysis command.

Use this when:

- you need a MAGMA step 2 option that is not wrapped as a named R argument
- you want to pass through an uncommon manual option directly

## How To Think About the Four Germline Argument Groups

`prepare_germline_scores()` exposes four MAGMA argument bundles:

- `gene_step1_args`
- `gene_step2_args`
- `reg_step1_args`
- `reg_step2_args`

These can be understood in either of two ways.

### Grouped by Target

- gene arguments = `gene_step1_args` + `gene_step2_args`
- regulatory arguments = `reg_step1_args` + `reg_step2_args`

### Grouped by MAGMA Stage

- step 1 arguments = `gene_step1_args` + `reg_step1_args`
- step 2 arguments = `gene_step2_args` + `reg_step2_args`

Both views are correct.

## What Is Not Yet Fully Wrapped

The package currently does **not** provide first-class named R parameters for
every MAGMA manual option.

If a user needs a MAGMA option that is not listed above, the intended path is:

1. pass it through `extra_args`
2. if it proves commonly useful, promote it later to a named wrapper argument

This keeps the API manageable while still preserving access to the broader
MAGMA CLI.
