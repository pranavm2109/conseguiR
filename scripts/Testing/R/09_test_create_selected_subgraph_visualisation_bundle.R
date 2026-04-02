#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
  library(igraph)
})

source("scripts/Internals/R/09_create_selected_subgraph_visualisation_bundle.R")

default_selected_nodes_path <- "data/processed/gene_gene_selected_subgraph_nodes.tsv"
default_selected_edges_path <- "data/processed/gene_gene_selected_subgraph_edges.tsv"
default_selected_summary_path <- "data/processed/gene_gene_selected_subgraph_summary.tsv"
default_selected_plot_path <- "data/processed/gene_gene_selected_subgraph_plot.pdf"
default_selected_bundle_test_output_dir <- "data/processed/test_outputs/selected_subgraph_bundle"

make_selected_bundle_test_path <- function(stem, ext = "") {
  dir.create(default_selected_bundle_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_selected_bundle_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

make_live_selected_subgraph_fixture <- function() {
  list(
    nodes = read_selected_subgraph_nodes(default_selected_nodes_path),
    edges = read_selected_subgraph_edges(default_selected_edges_path),
    summary = read_selected_subgraph_summary(default_selected_summary_path)
  )
}

test_prepare_selected_subgraph_visualisation_bundle_live_outputs <- function(print_bundle = TRUE) {
  bundle <- prepare_selected_subgraph_visualisation_bundle(
    nodes_path = default_selected_nodes_path,
    edges_path = default_selected_edges_path,
    summary_path = default_selected_summary_path
  )

  expect_true(is.data.frame(bundle$nodes))
  expect_true(is.data.frame(bundle$edges))
  expect_true(is.data.frame(bundle$summary))
  expect_true(inherits(bundle$igraph, "igraph"))
  expect_true(all(c("label", "node_class") %in% names(bundle$nodes)))
  expect_true(all(c("from", "to") %in% names(bundle$edges)))
  expect_true(all(c("id", "shared.name", "name") %in% names(bundle$cytoscape_nodes)))
  expect_true(all(c("source", "target", "interaction") %in% names(bundle$cytoscape_edges)))
  expect_equal(vcount(bundle$igraph), nrow(bundle$nodes))
  expect_equal(ecount(bundle$igraph), nrow(bundle$edges))
  expect_true(all(is.finite(bundle$nodes$prize)))

  if (isTRUE(print_bundle)) {
    message("Top bundled subgraph nodes for plotting:")
    print(bundle$nodes[order(-prize)][1:min(10L, .N), .(node_id, gene_name, prize, post_norm)])
  }

  invisible(bundle)
}

test_save_selected_subgraph_visualisation_bundle_outputs <- function() {
  bundle <- prepare_selected_subgraph_visualisation_bundle(
    nodes_path = default_selected_nodes_path,
    edges_path = default_selected_edges_path,
    summary_path = default_selected_summary_path
  )
  output_prefix <- make_selected_bundle_test_path("selected_subgraph_plot_bundle")

  save_selected_subgraph_visualisation_bundle(
    bundle = bundle,
    output_prefix = output_prefix,
    save_rds = TRUE,
    save_tables = TRUE
  )

  expect_true(file.exists(paste0(output_prefix, ".rds")))
  expect_true(file.exists(paste0(output_prefix, "_nodes.tsv.gz")))
  expect_true(file.exists(paste0(output_prefix, "_edges.tsv.gz")))
  expect_true(file.exists(paste0(output_prefix, "_summary.tsv")))
  expect_true(file.exists(paste0(output_prefix, "_cytoscape_nodes.tsv.gz")))
  expect_true(file.exists(paste0(output_prefix, "_cytoscape_edges.tsv.gz")))
}

test_save_selected_subgraph_plot_on_real_outputs <- function(
  file_path = default_selected_plot_path
) {
  bundle <- prepare_selected_subgraph_visualisation_bundle(
    nodes_path = default_selected_nodes_path,
    edges_path = default_selected_edges_path,
    summary_path = default_selected_summary_path
  )

  save_selected_subgraph_plot(
    bundle = bundle,
    file_path = file_path,
    width = 12,
    height = 10,
    dpi = 300
  )

  expect_true(file.exists(file_path))
  expect_gt(file.info(file_path)$size, 0)
  message("Saved selected subgraph plot to: ", file_path)
}

test_read_selected_subgraph_negative_missing_file <- function() {
  expect_error(
    read_selected_subgraph_nodes("data/processed/does_not_exist_nodes.tsv"),
    expected = "Selected subgraph node file does not exist",
    fixed = TRUE
  )

  expect_error(
    read_selected_subgraph_edges("data/processed/does_not_exist_edges.tsv"),
    expected = "Selected subgraph edge file does not exist",
    fixed = TRUE
  )

  expect_error(
    read_selected_subgraph_summary("data/processed/does_not_exist_summary.tsv"),
    expected = "Selected subgraph summary file does not exist",
    fixed = TRUE
  )
}

test_validate_selected_subgraph_negative_bad_nodes <- function() {
  expect_error(
    validate_selected_subgraph_nodes(data.table(node_id = "MYC", prize = 1)),
    expected = "Selected subgraph node table is missing required columns",
    fixed = TRUE
  )

  expect_error(
    validate_selected_subgraph_nodes(data.table(
      node_id = c("MYC", "MYC"),
      gene_name = c("MYC", "MYC"),
      prize = c(1, 2)
    )),
    expected = "duplicated `node_id` values",
    fixed = TRUE
  )

  expect_error(
    validate_selected_subgraph_nodes(data.table(
      node_id = c("MYC", "BCL6"),
      gene_name = c("MYC", "BCL6"),
      prize = c("bad", "1.2")
    )),
    expected = "non-numeric `prize` values",
    fixed = TRUE
  )
}

test_validate_selected_subgraph_negative_bad_edges <- function() {
  expect_error(
    validate_selected_subgraph_edges(data.table(gene_u = "MYC")),
    expected = "Selected subgraph edge table is missing required columns",
    fixed = TRUE
  )

  expect_error(
    validate_selected_subgraph_edges(data.table(gene_u = "MYC", gene_v = "")),
    expected = "contains missing gene endpoints",
    fixed = TRUE
  )
}

test_validate_selected_subgraph_negative_bad_summary <- function() {
  expect_error(
    validate_selected_subgraph_summary(data.table(target_genes = 50)),
    expected = "Selected subgraph summary table is missing required columns",
    fixed = TRUE
  )
}

test_prepare_subgraph_plot_edges_negative_no_overlap <- function() {
  fixture <- make_live_selected_subgraph_fixture()
  bad_edges <- data.table(gene_u = "NOT_IN_GRAPH_1", gene_v = "NOT_IN_GRAPH_2")

  expect_error(
    prepare_subgraph_plot_edges(bad_edges, fixture$nodes),
    expected = "do not overlap the selected node identifiers",
    fixed = TRUE
  )
}

test_save_selected_subgraph_plot_negative_missing_file_path <- function() {
  bundle <- prepare_selected_subgraph_visualisation_bundle(
    nodes_path = default_selected_nodes_path,
    edges_path = default_selected_edges_path,
    summary_path = default_selected_summary_path
  )

  expect_error(
    save_selected_subgraph_plot(bundle = bundle, file_path = ""),
    expected = "`file_path` must be provided",
    fixed = TRUE
  )
}

main <- function() {
  test_that("selected subgraph plot bundle builds from the real end-to-end outputs", {
    test_prepare_selected_subgraph_visualisation_bundle_live_outputs(print_bundle = TRUE)
  })

  test_that("selected subgraph plot bundle outputs can be saved", {
    test_save_selected_subgraph_visualisation_bundle_outputs()
  })

  test_that("selected subgraph plot can be saved from the real end-to-end outputs", {
    test_save_selected_subgraph_plot_on_real_outputs()
  })

  test_that("selected subgraph readers fail clearly for missing files", {
    test_read_selected_subgraph_negative_missing_file()
  })

  test_that("selected subgraph node validators fail clearly for malformed node tables", {
    test_validate_selected_subgraph_negative_bad_nodes()
  })

  test_that("selected subgraph edge validators fail clearly for malformed edge tables", {
    test_validate_selected_subgraph_negative_bad_edges()
  })

  test_that("selected subgraph summary validators fail clearly for malformed summaries", {
    test_validate_selected_subgraph_negative_bad_summary()
  })

  test_that("plot-edge preparation fails clearly when edges do not overlap nodes", {
    test_prepare_subgraph_plot_edges_negative_no_overlap()
  })

  test_that("plot saving fails clearly when file_path is missing", {
    test_save_selected_subgraph_plot_negative_missing_file_path()
  })
}

if (sys.nframe() == 0) {
  main()
}
