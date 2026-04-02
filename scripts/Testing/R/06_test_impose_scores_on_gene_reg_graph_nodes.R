#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
  library(igraph)
})

source("scripts/Internals/R/06_impose_scores_on_gene_reg_graph_nodes.R")

default_graph_rds_path <- "data/processed/gene_reg_graph_no_scores.rds"
default_gene_reg_scored_test_output_dir <- "data/processed/test_outputs/gene_reg_scored"

make_gene_reg_scored_test_path <- function(stem, ext = "") {
  dir.create(default_gene_reg_scored_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_gene_reg_scored_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

make_live_gene_reg_score_fixture <- function() {
  graph <- read_gene_reg_graph_no_scores(default_graph_rds_path)
  nodes <- extract_gene_reg_graph_nodes(graph)

  gene_ids <- nodes[node_type == "gene", unique(node_id)][1:3]
  reg_ids <- nodes[node_type == "reg", unique(node_id)][1:3]

  list(
    graph = graph,
    nodes = nodes,
    gene_somatic_scores = data.table(
      gene_id = gene_ids,
      zstat = c(2.1, -1.4, 0.8)
    ),
    gene_germline_scores = data.table(
      gene_id = gene_ids,
      zstat = c(1.2, 0.4, -0.6)
    ),
    reg_somatic_scores = data.table(
      reg_elem_id = reg_ids,
      zstat = c(0.9, -0.3, 1.5)
    ),
    reg_germline_scores = data.table(
      reg_elem_id = reg_ids,
      zstat = c(-0.5, 0.7, 1.1)
    ),
    reg_epigenomic_scores = data.table(
      reg_elem_id = reg_ids,
      zstat = c(3.0, -2.0, 0.5)
    )
  )
}

test_impose_scores_on_gene_reg_nodes_live_fixture <- function() {
  fixture <- make_live_gene_reg_score_fixture()

  scored_nodes <- impose_scores_on_gene_reg_nodes(
    nodes = fixture$nodes,
    gene_somatic_scores = fixture$gene_somatic_scores,
    gene_germline_scores = fixture$gene_germline_scores,
    reg_somatic_scores = fixture$reg_somatic_scores,
    reg_germline_scores = fixture$reg_germline_scores,
    reg_epigenomic_scores = fixture$reg_epigenomic_scores
  )

  expect_true(is.data.frame(scored_nodes))
  expect_true(all(c("somatic_score", "germline_score", "epigenomic_score") %in% names(scored_nodes)))

  gene_nodes <- scored_nodes[node_type == "gene"]
  reg_nodes <- scored_nodes[node_type == "reg"]

  expect_true(all(gene_nodes$epigenomic_score == 0))
  expect_true(any(reg_nodes$epigenomic_score != 0))
  expect_true(any(gene_nodes$somatic_score != 0))
  expect_true(any(gene_nodes$germline_score != 0))
  invisible(scored_nodes)
}

test_impose_scores_on_gene_reg_graph_live_fixture <- function() {
  fixture <- make_live_gene_reg_score_fixture()

  result <- impose_scores_on_gene_reg_graph(
    graph = fixture$graph,
    gene_somatic_scores = fixture$gene_somatic_scores,
    gene_germline_scores = fixture$gene_germline_scores,
    reg_somatic_scores = fixture$reg_somatic_scores,
    reg_germline_scores = fixture$reg_germline_scores,
    reg_epigenomic_scores = fixture$reg_epigenomic_scores
  )

  expect_true(inherits(result$graph, "igraph"))
  expect_true(is.data.frame(result$nodes))
  expect_true(is.data.frame(result$edges))
  expect_true(all(c("somatic_score", "germline_score", "epigenomic_score") %in% vertex_attr_names(result$graph)))
  invisible(result)
}

test_prepare_scored_gene_reg_graph_saves_outputs <- function(print_scores = TRUE) {
  fixture <- make_live_gene_reg_score_fixture()
  output_prefix <- make_gene_reg_scored_test_path("gene_reg_graph_scored_test")

  result <- prepare_scored_gene_reg_graph(
    graph = fixture$graph,
    output_prefix = output_prefix,
    gene_somatic_scores = fixture$gene_somatic_scores,
    gene_germline_scores = fixture$gene_germline_scores,
    reg_somatic_scores = fixture$reg_somatic_scores,
    reg_germline_scores = fixture$reg_germline_scores,
    reg_epigenomic_scores = fixture$reg_epigenomic_scores,
    save_outputs = TRUE
  )

  expect_true(file.exists(paste0(output_prefix, ".rds")))
  expect_true(file.exists(paste0(output_prefix, "_nodes.tsv.gz")))
  expect_true(file.exists(paste0(output_prefix, "_edges.tsv.gz")))

  if (isTRUE(print_scores)) {
    message("Scored gene-reg nodes with injected non-zero scores:")
    print(
      result$nodes[
        node_id %in% c(
          fixture$gene_somatic_scores$gene_id,
          fixture$reg_somatic_scores$reg_elem_id
        ),
        .(node_id, node_type, somatic_score, germline_score, epigenomic_score)
      ]
    )
  }

  invisible(result)
}

test_standardize_score_tables_negative_cases <- function() {
  expect_error(
    standardize_gene_score_table(data.table(gene = "TP53", z = 1), "somatic"),
    expected = "Gene somatic score table is missing required columns",
    fixed = TRUE
  )

  expect_error(
    standardize_reg_score_table(data.table(reg = "GH01J000013", z = 1), "epigenomic"),
    expected = "Regulatory-element epigenomic score table is missing required columns",
    fixed = TRUE
  )
}

test_validate_gene_reg_graph_nodes_negative_cases <- function() {
  expect_error(
    validate_gene_reg_graph_nodes(data.table(node_id = "TP53", node_type = "gene")),
    expected = "Gene-reg node table is missing required columns",
    fixed = TRUE
  )

  expect_error(
    validate_gene_reg_graph_nodes(
      data.table(
        name = "TP53",
        node_id = "TP53",
        node_type = "protein"
      )
    ),
    expected = "unsupported `node_type` values",
    fixed = TRUE
  )
}

test_read_gene_reg_graph_no_scores_negative_missing_file <- function() {
  expect_error(
    read_gene_reg_graph_no_scores("data/processed/does_not_exist.rds"),
    expected = "Gene-reg graph file does not exist",
    fixed = TRUE
  )
}

test_prepare_scored_gene_reg_graph_negative_bad_gene_scores <- function() {
  fixture <- make_live_gene_reg_score_fixture()

  expect_error(
    prepare_scored_gene_reg_graph(
      graph = fixture$graph,
      gene_somatic_scores = data.table(gene = "TP53", z = 1),
      save_outputs = FALSE
    ),
    expected = "Gene somatic score table is missing required columns",
    fixed = TRUE
  )
}

test_prepare_scored_gene_reg_graph_negative_bad_reg_scores <- function() {
  fixture <- make_live_gene_reg_score_fixture()

  expect_error(
    prepare_scored_gene_reg_graph(
      graph = fixture$graph,
      reg_epigenomic_scores = data.table(reg = "GH01J000013", z = 1),
      save_outputs = FALSE
    ),
    expected = "Regulatory-element epigenomic score table is missing required columns",
    fixed = TRUE
  )
}

main <- function() {
  test_that("gene-reg node scoring works on the live backend node table", {
    test_impose_scores_on_gene_reg_nodes_live_fixture()
  })

  test_that("gene-reg graph scoring updates igraph vertex attributes", {
    test_impose_scores_on_gene_reg_graph_live_fixture()
  })

  test_that("scored gene-reg graph outputs can be saved", {
    test_prepare_scored_gene_reg_graph_saves_outputs(print_scores = TRUE)
  })

  test_that("score-table validators report clear errors for malformed score tables", {
    test_standardize_score_tables_negative_cases()
  })

  test_that("gene-reg node validators report clear errors for broken node tables", {
    test_validate_gene_reg_graph_nodes_negative_cases()
  })

  test_that("reading the backend no-score graph fails clearly when the file is missing", {
    test_read_gene_reg_graph_no_scores_negative_missing_file()
  })

  test_that("scored graph preparation fails clearly for malformed gene score tables", {
    test_prepare_scored_gene_reg_graph_negative_bad_gene_scores()
  })

  test_that("scored graph preparation fails clearly for malformed regulatory score tables", {
    test_prepare_scored_gene_reg_graph_negative_bad_reg_scores()
  })
}

if (sys.nframe() == 0) {
  main()
}
