# Bioconductor Readiness Audit

This document tracks the package against common Bioconductor expectations from
the official contribution guide.

## Already improved

- Exported functions have man pages and runnable examples.
- A narrative vignette exists in `vignettes/conseguiR-overview.Rmd`.
- A second report-style analysis document exists in
  `vignettes/conseguiR-analysis-report.Rmd`.
- Package startup is now lightweight and no longer auto-initializes backend
  graphs or auto-discovers Python runtimes.
- Python-backed stages now run through a managed `basilisk` environment rather
  than a hard-coded user conda environment.
- Standard package files `NEWS.md` and `INSTALL` are now present.
- Generated vignette markdown files are excluded from the build.
- A standard `tests/testthat` entry point is now present.

## Remaining likely blockers

### 1. Python runtime strategy

The diffusion and selected-subgraph stages still rely on Python, but they now
use a managed `basilisk` environment instead of a user-managed interpreter.
This substantially improves Bioconductor alignment, but it still needs a clean
end-to-end validation pass on a fresh machine / build context.

Recommended fix:

- validate the `basilisk` environment provisioning path on a clean setup
- confirm the managed environment works under package build/check conditions

### 2. Test layout and scope

The historical tests still live under `scripts/Testing/R` and
`scripts/Testing/Python`. A minimal `tests/testthat` harness now exists, but
the full package test suite still needs to be migrated into standard package
testing layout.

Recommended fix:

- move or rewrite the existing script-style tests under `tests/testthat/`
- make sure package tests do not depend on local development-only outputs under
  `data/processed/`

### 3. Package data footprint

The repository contains large development data outside the build, and
`inst/extdata` is still sizeable. This should be reviewed before submission.

Recommended fix:

- keep only the smallest necessary example files in `inst/extdata`
- move large supporting resources to an external download/cache path or a hub
  strategy if they are essential

### 4. Metadata polish

The package now has `biocViews`, but it still likely needs:

- final submission-ready versioning
- project `URL` and `BugReports`
- a `CITATION` file if the package/paper citation should be explicit

## Practical next steps

1. migrate the script-style tests into `tests/testthat`
2. decide whether Python-backed stages are optional or `basilisk`-managed
3. trim `inst/extdata` to submission-friendly example data
4. add final package metadata (`URL`, `BugReports`, `CITATION`)
5. run `R CMD check` and compare against the Bioconductor build expectations
