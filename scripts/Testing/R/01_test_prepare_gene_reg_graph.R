#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(testthat)
})

source("scripts/Internals/R/01_prepare_gene_reg_graph.R")

genehancer_interactions_path <- "data/raw/GeneHancer/gh_interactions_hg38_primary_assembly"
genehancer_reg_elements_path <- "data/raw/GeneHancer/gh_reg_elements_hg38_primary_assembly"
default_gene_reg_test_output_dir <- "data/processed/test_outputs/gene_reg_graph"

make_gene_reg_test_path <- function(stem, ext = "") {
  dir.create(default_gene_reg_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_gene_reg_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

make_gene_reg_test_inputs <- function() {
  interactions <- fread(genehancer_interactions_path, nrows = 3000)
  reg_ids <- unique(interactions$geneHancerIdentifier)
  reg_elements <- fread(genehancer_reg_elements_path)
  reg_elements <- reg_elements[name %in% reg_ids]

  interactions_file <- make_gene_reg_test_path("interactions", ".tsv")
  reg_elements_file <- make_gene_reg_test_path("reg_elements", ".tsv")

  fwrite(interactions, interactions_file, sep = "\t")
  fwrite(reg_elements, reg_elements_file, sep = "\t")

  list(
    interactions_path = interactions_file,
    reg_elements_path = reg_elements_file
  )
}

test_prepare_gene_reg_graph <- function() {
  inputs <- make_gene_reg_test_inputs()
  output_prefix <- make_gene_reg_test_path("gene_reg_graph_test")

  result <- prepare_gene_reg_graph(config = list(
    interactions_path = inputs$interactions_path,
    reg_elements_path = inputs$reg_elements_path,
    output_prefix = output_prefix,
    min_link_value = 0,
    keep_self_loops = FALSE,
    directed = FALSE
  ))

  expect_s3_class(result$graph, "igraph")
  expect_gt(nrow(result$nodes), 0L)
  expect_gt(nrow(result$edges), 0L)
  expect_true(all(c("name", "node_id", "node_type") %in% names(result$nodes)))
  expect_true(all(c("from", "to", "weight", "confidence") %in% names(result$edges)))
  expect_true(all(result$nodes$node_type %in% c("gene", "reg")))

  node_type_map <- setNames(result$nodes$node_type, result$nodes$name)
  from_types <- unname(node_type_map[result$edges$from])
  to_types <- unname(node_type_map[result$edges$to])

  expect_true(all(from_types == "reg"))
  expect_true(all(to_types == "gene"))

  expected_files <- c(
    paste0(output_prefix, ".rds"),
    paste0(output_prefix, "_nodes.tsv.gz"),
    paste0(output_prefix, "_edges.tsv.gz")
  )
  expect_true(all(file.exists(expected_files)))

  message("Gene-reg graph test passed.")
  invisible(result)
}

main <- function() {
  test_that("gene-reg graph preparation works", {
    result <- test_prepare_gene_reg_graph()
    expect_s3_class(result$graph, "igraph")
  })
}

if (sys.nframe() == 0) {
  main()
}
