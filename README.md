# conseguiR

`conseguiR` is an R package skeleton for convergent germline, somatic, and regulatory integration on graphs.

## Included structure

- `R/` exported and internal function stubs
- `inst/extdata/STRING/` placeholder location for STRING resources
- `inst/extdata/GeneHancer/` placeholder location for GeneHancer resources
- `inst/extdata/example_data/` placeholder location for example inputs
- `data-raw/` scripts/placeholders for backend preparation
- `tests/testthat/` package test skeleton

## Notes

This scaffold is meant to be filled in with real implementations. Functions are documented and organized into:

- user-facing exported API
- internal helper functions
- pipeline wrappers

## Suggested next steps

1. Run `devtools::document()` to generate `NAMESPACE` and `.Rd` files.
2. Add backend resources under `inst/extdata/`.
3. Implement scoring wrappers around dndscv, fishHook, and MAGMA.
4. Add tests and example datasets.
