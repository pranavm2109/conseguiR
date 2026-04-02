#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
})

source("scripts/Internals/R/08_call_subgraph.R")
default_step8_test_output_dir <- "data/processed/test_outputs/step8"

make_step8_test_path <- function(stem, ext = "") {
  dir.create(default_step8_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_step8_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

python_runtime_candidates <- function() {
  candidates <- c("python", "python3")
  candidates[nzchar(Sys.which(candidates))]
}

python_module_available <- function(module_name) {
  candidates <- python_runtime_candidates()
  if (length(candidates) == 0L) {
    return(FALSE)
  }

  for (python_cmd in candidates) {
    status <- system2(
      python_cmd,
      args = c("-c", shQuote(paste0("import ", module_name))),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_status <- attr(status, "status")
    if (is.null(exit_status) || identical(exit_status, 0L)) {
      return(TRUE)
    }
  }
  FALSE
}

make_diffusion_fixture <- function() {
  data.table(
    node_id = c("G1", "G2", "G3"),
    gene_name = c("G1", "G2", "G3"),
    post_norm = c(1.0, 1.5, 0.8),
    post_germline = c(0.8, 1.0, 0.5),
    post_somatic = c(0.2, 0.3, 0.2),
    post_epigenomic = c(0.5, 0.6, 0.4)
  )
}

make_gene_gene_graph_fixture <- function() {
  list(
    nodes = data.table(
      node_id = c("G1", "G2", "G3"),
      node_type = c("gene", "gene", "gene")
    ),
    edges = data.table(
      from = c("G1", "G2", "G2"),
      to = c("G2", "G3", "G1"),
      confidence = c(0.9, 0.8, 0.7),
      weight = c(1.0, 1.0, 1.0)
    )
  )
}

python_ortools_available <- function() {
  candidates <- python_runtime_candidates()
  if (length(candidates) == 0L) {
    return(FALSE)
  }

  for (python_cmd in candidates) {
    status <- system2(
      python_cmd,
      args = c("-c", shQuote("import ortools")),
      stdout = TRUE,
      stderr = TRUE
    )
    exit_status <- attr(status, "status")
    if (is.null(exit_status) || identical(exit_status, 0L)) {
      return(TRUE)
    }
  }
  FALSE
}

test_validate_diffusion_results <- function() {
  diffusion <- make_diffusion_fixture()
  validated <- validate_diffusion_results(diffusion, prize_column = "post_norm")

  expect_true(is.data.table(validated))
  expect_equal(nrow(validated), 3L)
}

test_validate_diffusion_results_negative <- function() {
  bad_diffusion <- data.table(node_id = "G1", post_norm = 1.0)
  expect_error(
    validate_diffusion_results(bad_diffusion, prize_column = "post_norm"),
    expected = "Diffusion results are missing required columns",
    fixed = TRUE
  )
}

test_validate_gene_gene_graph_inputs <- function() {
  fixture <- make_gene_gene_graph_fixture()
  validated_nodes <- validate_gene_gene_nodes(fixture$nodes)
  validated_edges <- validate_gene_gene_edges(fixture$edges)

  expect_true(is.data.table(validated_nodes))
  expect_true(is.data.table(validated_edges))
}

test_validate_gene_gene_graph_inputs_negative <- function() {
  bad_nodes <- data.table(node_id = "G1", node_type = "protein")
  expect_error(
    validate_gene_gene_nodes(bad_nodes),
    expected = "Gene-gene node table contains unsupported node_type values",
    fixed = TRUE
  )

  bad_edges <- data.table(from = "G1", to = "G2", weight = 1.0)
  expect_error(
    validate_gene_gene_edges(bad_edges),
    expected = "Gene-gene edge table is missing required columns",
    fixed = TRUE
  )
}

test_run_cardinality_subgraph_calling_integration <- function() {
  if (length(python_runtime_candidates()) == 0L) {
    skip("Python is not available on PATH.")
  }
  if (!python_module_available("numpy") || !python_module_available("pandas")) {
    skip("Required Python packages are not installed for step 8 integration testing.")
  }
  if (!python_ortools_available()) {
    skip("OR-Tools is not available for live step 8 testing.")
  }

  output_dir <- make_step8_test_path("conseguiR_step8")
  result <- run_cardinality_subgraph_calling(
    output_dir = output_dir,
    output_stem = "gene_gene_selected_subgraph_test",
    target_genes = 2L,
    candidate_pool_size = 3L,
    max_edges_in_model = 50L,
    max_time_seconds = 30L,
    num_workers = 1L
  )

  expect_true(file.exists(result$output_paths$nodes_path))
  expect_true(file.exists(result$output_paths$edges_path))
  expect_true(file.exists(result$output_paths$summary_path))
  expect_true(file.exists(result$output_paths$graphml_path))
}

main <- function() {
  test_that("diffusion result validation works", {
    test_validate_diffusion_results()
  })

  test_that("diffusion validation reports missing required columns", {
    test_validate_diffusion_results_negative()
  })

  test_that("gene-gene graph validation works", {
    test_validate_gene_gene_graph_inputs()
  })

  test_that("gene-gene graph validation rejects malformed input", {
    test_validate_gene_gene_graph_inputs_negative()
  })

  test_that("step 8 Python-backed subgraph wrapper runs when OR-Tools is available", {
    test_run_cardinality_subgraph_calling_integration()
  })
}

if (sys.nframe() == 0) {
  main()
}
