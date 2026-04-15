# Bioconductor Submission Checklist

This document tracks `conseguiR` against the Bioconductor submissions guidance,
especially:

- Chapter 1 "Bioconductor Package Submissions"
- the linked package-development chapters on `README`, `DESCRIPTION`,
  `NAMESPACE`, `LICENSE`, `INSTALL`, documentation, package data, unit tests,
  R code, third-party code, and `.gitignore`

Primary reference:

- https://contributions.bioconductor.org/bioconductor-package-submissions.html

The goal of this file is practical triage:

- what is already in good shape
- what is partially addressed but still review-sensitive
- what still needs to be fixed before submission

## Done

These items are either already aligned with the submission page or have been
brought into much better shape during the recent cleanup.

### Scope and package type

- `conseguiR` is clearly a software package aimed at genomic / regulatory /
  mutation-analysis workflows, so it fits the broad Bioconductor scope for
  high-throughput biological data analysis.
- `biocViews` is present in `DESCRIPTION`.

### Core package metadata and standard files

- `README.md` exists.
- `DESCRIPTION` exists and now passes the `DESCRIPTION meta-information` stage
  in staged `R CMD check`.
- `NEWS.md` exists.
- `INSTALL` exists.
- `LICENSE` now uses a valid MIT DCF stub:
  - `YEAR: 2026`
  - `COPYRIGHT HOLDER: Pranav Mishra`

### Documentation baseline

- Exported user-facing functions have man pages.
- A narrative vignette exists in `vignettes/conseguiR-overview.Rmd`.
- The package documentation is no longer failing at the `DESCRIPTION`
  meta-information stage because of a malformed `LICENSE`.

### Runtime and package behavior

- Package startup is lightweight and no longer tries to initialize heavy
  backend state automatically.
- Python-backed stages now use a managed `basilisk` environment instead of a
  hard-coded user conda setup.
- The bundled MAGMA executable has been removed from the built package tarball,
  so the old "undeclared executable file" warning is fixed.
- MAGMA is now treated as an external prerequisite, with resolution via:
  - explicit `magma_path`
  - `options(conseguiR.magma_path = ...)`
  - `CONSEGUIR_MAGMA_PATH`
  - `magma` on `PATH`
  - a repository-local development fallback when present

### Recent check issues already fixed

- The `setNames` "no visible global function definition" note was fixed by
  using `stats::setNames()`.
- The earlier "unused Imports" note was resolved in staged checking by adding
  targeted `@importFrom` declarations.
- The repository-local `man/initialize_backend_graphs.Rd` now documents all of
  the function arguments.
- macOS `._*` shadow files were removed from `R/`, which had been interfering
  with roxygen behavior.

## Partially Done

These areas are moving in the right direction, but they are still likely to
attract reviewer attention or still need one more round of hardening.

### Interoperability with Bioconductor infrastructure

- The package uses Bioconductor-relevant packages and genomic concepts.
- However, the user-facing API is still mainly built around package-specific
  bundles and sourced script layers rather than deeper integration with
  canonical Bioconductor container classes throughout the workflow.
- This does not automatically disqualify the package, but reviewers may ask
  for a clearer interoperability story.

### Documentation quality

- The package has a vignette and function-level documentation.
- Documentation still needs one more pass for submission polish:
  - examples of expected input formats
  - clearer explanation of external prerequisites
  - a cleaner narrative of how the package interoperates with Bioconductor
    infrastructure

### Python / third-party code story

- The move to `basilisk` is a major improvement.
- The package still has meaningful Python-backed execution and external MAGMA
  dependence, so the "third-party code" review area is still sensitive.
- What remains is less about architecture and more about proving that the
  resulting user experience is robust in a clean build/check environment.

### Testing

- A `tests/testthat` harness exists.
- Historical tests still largely live under `scripts/Testing`.
- This means the package has tests, but not yet in a submission-polished,
  standard package layout.

### Build/check progress

- Staged `R CMD build` and `R CMD check` are much healthier than before.
- The package now passes the following previously failing check stages in staged
  runs:
  - `checking for executable files`
  - `checking DESCRIPTION meta-information`
  - `checking dependencies in R code`
- We still have not yet driven a full submission-style run to a genuinely clean
  end state under devel Bioconductor plus `BiocCheck(new-package = TRUE)`.

## Must Fix Before Submission

These are the items that should be treated as actual submission blockers until
resolved.

### 1. Clean default branch / package-only repository state

The submissions page says the default branch used for submission must contain
only package code. In practice, `conseguiR` still looks like an active research
repository rather than a package-only submission branch.

Needs work:

- remove or relocate development-only artifacts from the submission branch
- ensure the default branch is package-focused
- verify that no extra non-package files bleed into `R CMD build`

### 2. Full `R CMD check` and `BiocCheck` on the right environment

Bioconductor explicitly says the package will be subjected to:

- `R CMD build`
- `R CMD check`
- `BiocCheckGitClone()`
- `BiocCheck(new-package = TRUE)`

Needs work:

- run `BiocCheck(new-package = TRUE)` and fix all actionable issues
- run the package against the devel version of Bioconductor
- verify behavior on Linux, macOS, and Windows expectations as far as possible

### 3. Final package metadata polish

Still missing or incomplete:

- `URL` field in `DESCRIPTION`
- `BugReports` field in `DESCRIPTION`
- `CITATION` file
- final submission-ready versioning / release discipline

### 4. Test migration into standard package layout

Bioconductor expects a conventional, maintained testing story.

Needs work:

- migrate or rewrite meaningful tests under `tests/testthat`
- reduce reliance on `scripts/Testing`
- ensure tests are self-contained and do not depend on local development output

### 5. Package data footprint

`inst/extdata` is still large enough to merit attention.

Current staged-check signal:

- installed size about `43.1 MB`
- `extdata` about `42.5 MB`

Needs work:

- decide which resources are truly package-owned examples
- move anything larger or more reference-like to a more appropriate strategy if
  possible
- document why each shipped resource needs to live inside the package

### 6. Third-party dependency clarity

Even though MAGMA externalization is fixed technically, the submission story
still needs to be very clear.

Needs work:

- document MAGMA as an optional / external prerequisite as cleanly as possible
- ensure package functions fail gracefully when MAGMA is absent
- make sure the package is still useful and checkable even when the external
  tool is not installed

## Current High-Priority Queue

If we were prioritizing strictly for submission readiness, the next steps
should be:

1. Add `URL`, `BugReports`, and a `CITATION` file.
2. Run `BiocCheck(new-package = TRUE)` and capture every issue.
3. Clean the default branch so it looks like package source, not a mixed
   package-plus-development workspace.
4. Migrate the important tests into `tests/testthat`.
5. Reassess `inst/extdata` and shrink it where possible.
6. Do one final pass on Python / MAGMA dependency messaging and failure modes.

## Notes on Recent Check Results

Interpreting recent staged checks requires care because a few notes were caused
by temporary staging mistakes rather than the repository itself.

### Real issues that were fixed

- bundled MAGMA executable warning
- invalid MIT license stub
- `setNames` namespace note
- broad `Imports` note

### Temporary staging artifacts, not current repo blockers

- "left-over files" notes caused by copying `man/*.Rd` into the staged root
  instead of into `man/`
- stale staged copies of `initialize_backend_graphs.Rd`

## Definition of "Submission Ready"

For this repository, we should consider `conseguiR` submission-ready only when
all of the following are true:

- default branch contains only package-relevant source
- staged `R CMD build` and `R CMD check` are clean aside from unavoidable
  environment-only limitations
- `BiocCheck(new-package = TRUE)` is clean or all remaining issues have
  deliberate, documented justification
- metadata is complete (`URL`, `BugReports`, `CITATION`)
- tests are in standard package layout
- shipped data and third-party dependencies have a defensible Bioconductor
  story
