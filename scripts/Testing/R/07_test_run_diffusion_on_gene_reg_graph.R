#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
  devtools::load_all(".", quiet = TRUE)
})

source("scripts/Internals/R/06b_python_basilisk.R")
source("scripts/Internals/R/07_run_diffusion_on_gene_reg_graph.R")
default_step7_test_output_dir <- "data/processed/test_outputs/step7"

make_step7_test_path <- function(stem, ext = "") {
  dir.create(default_step7_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_step7_test_output_dir,
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
generate_scored_gene_reg_graph_fixtures <- function() {
  nodes <- data.table(
    name = c("G1", "G2", "R1", "R2"),
    node_id = c("G1", "G2", "R1", "R2"),
    node_type = c("gene", "gene", "reg", "reg"),
    somatic_score = c(1.0, -0.5, 0.2, 0.3),
    germline_score = c(0.5, 0.0, 0.1, -0.2),
    epigenomic_score = c(0.0, 0.0, 1.5, -0.8)
  )

  edges <- data.table(
    from = c("R1", "G1", "R2", "G2"),
    to = c("G1", "R1", "G2", "R2"),
    confidence = c(0.9, 0.8, 0.7, 0.85)
  )

  list(nodes = nodes, edges = edges)
}

write_scored_gene_reg_fixture_files <- function() {
  fixture <- generate_scored_gene_reg_graph_fixtures()
  output_prefix <- make_step7_test_path("scored_gene_reg_fixture")
  nodes_path <- paste0(output_prefix, "_nodes.tsv")
  edges_path <- paste0(output_prefix, "_edges.tsv")
  fwrite(fixture$nodes, nodes_path, sep = "\t")
  fwrite(fixture$edges, edges_path, sep = "\t")
  list(
    nodes_path = nodes_path,
    edges_path = edges_path
  )
}

test_validate_scored_gene_reg_nodes <- function() {
  fixture <- generate_scored_gene_reg_graph_fixtures()
  validated <- validate_scored_gene_reg_nodes(fixture$nodes)

  expect_true(is.data.table(validated))
  expect_equal(nrow(validated), 4L)
  expect_true(all(c("somatic_score", "germline_score", "epigenomic_score") %in% names(validated)))
}

test_validate_scored_gene_reg_nodes_negative_missing_col <- function() {
  bad_nodes <- data.table(
    node_id = c("G1", "R1"),
    node_type = c("gene", "reg"),
    somatic_score = c(1.0, 0.1)
  )

  expect_error(
    validate_scored_gene_reg_nodes(bad_nodes),
    regexp = "Scored gene-reg node table is missing required columns",
    fixed = TRUE
  )
}

test_validate_scored_gene_reg_nodes_negative_bad_type <- function() {
  bad_nodes <- data.table(
    node_id = c("G1", "X1"),
    node_type = c("gene", "protein"),
    somatic_score = c(1.0, 0.1),
    germline_score = c(0.5, 0.2),
    epigenomic_score = c(0.0, 0.0)
  )

  expect_error(
    validate_scored_gene_reg_nodes(bad_nodes),
    regexp = "Unsupported node_type values found in scored gene-reg nodes",
    fixed = TRUE
  )
}

test_validate_scored_gene_reg_edges_negative <- function() {
  bad_edges <- data.table(
    from = c("R1", NA),
    to = c("G1", "G2"),
    confidence = c(0.9, 0.8)
  )

  expect_error(
    validate_scored_gene_reg_edges(bad_edges),
    regexp = "Scored gene-reg edge table contains missing endpoint identifiers.",
    fixed = TRUE
  )
}

test_run_gene_reg_diffusion_integration <- function() {
  if (length(python_runtime_candidates()) == 0L) {
    skip("Python is not available on PATH.")
  }

  if (!python_module_available("numpy") || !python_module_available("pandas")) {
    skip("Required Python packages are not installed for step 7 integration testing.")
  }

  fixture_paths <- write_scored_gene_reg_fixture_files()
  output_dir <- make_step7_test_path("conseguiR_step7")
  result <- run_gene_reg_diffusion(
    nodes_path = fixture_paths$nodes_path,
    edges_path = fixture_paths$edges_path,
    output_dir = output_dir,
    output_stem = "gene_reg_graph_diffusion_test",
    top_n_to_save = 5L
  )

  expect_true(file.exists(result$output_paths$all_genes_path))
  expect_true(file.exists(result$output_paths$top_genes_path))
  expect_equal(result$config$output_dir, output_dir)
  diffusion_dt <- fread(result$output_paths$all_genes_path)
  expect_true(all(c("node_id", "gene_name", "post_integrated") %in% names(diffusion_dt)))
  expect_equal(nrow(diffusion_dt), 2L)
}

main <- function() {
  test_that("scored gene-reg node validation works", {
    test_validate_scored_gene_reg_nodes()
  })

  test_that("scored gene-reg node validation rejects bad input", {
    test_validate_scored_gene_reg_nodes_negative_missing_col()
    test_validate_scored_gene_reg_nodes_negative_bad_type()
  })

  test_that("scored gene-reg edge validation rejects missing endpoints", {
    test_validate_scored_gene_reg_edges_negative()
  })

  test_that("gate 7 Python-backed diffusion wrapper runs to completion", {
    test_run_gene_reg_diffusion_integration()
  })
}

if (sys.nframe() == 0) {
  main()
}
