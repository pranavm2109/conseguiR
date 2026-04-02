#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(testthat)
})

source("scripts/Internals/R/02_prepare_gene_gene_graph.R")

string_links_path <- "data/raw/STRING/9606.protein.links.v12.0.txt"
string_info_path <- "data/raw/STRING/9606.protein.info.v12.0.txt"
default_gene_gene_test_output_dir <- "data/processed/test_outputs/gene_gene_graph"

make_gene_gene_test_path <- function(stem, ext = "") {
  dir.create(default_gene_gene_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_gene_gene_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

make_gene_gene_test_inputs <- function() {
  links <- fread(string_links_path, nrows = 5000)
  protein_ids <- unique(c(links$protein1, links$protein2))
  protein_info <- fread(string_info_path)
  protein_info <- protein_info[`#string_protein_id` %in% protein_ids]

  links_file <- make_gene_gene_test_path("links", ".tsv")
  info_file <- make_gene_gene_test_path("info", ".tsv")

  fwrite(links, links_file, sep = "\t")
  fwrite(protein_info, info_file, sep = "\t")

  list(
    links_path = links_file,
    info_path = info_file
  )
}

test_prepare_gene_gene_graph <- function() {
  inputs <- make_gene_gene_test_inputs()
  output_prefix <- make_gene_gene_test_path("gene_gene_graph_test")

  result <- prepare_gene_gene_graph(config = list(
    protein_links_path = inputs$links_path,
    protein_info_path = inputs$info_path,
    output_prefix = output_prefix,
    min_combined_score = 400,
    directed = FALSE,
    collapse_to_gene_level = TRUE
  ))

  expect_s3_class(result$graph, "igraph")
  expect_gt(nrow(result$nodes), 0L)
  expect_gt(nrow(result$edges), 0L)
  expect_true(all(c("name", "node_id", "node_type") %in% names(result$nodes)))
  expect_true(all(c("from", "to", "confidence", "weight", "n_protein_edges") %in% names(result$edges)))
  expect_true(all(result$nodes$node_type == "gene"))
  expect_false(any(result$edges$from == result$edges$to))

  expected_files <- c(
    paste0(output_prefix, ".rds"),
    paste0(output_prefix, "_nodes.tsv.gz"),
    paste0(output_prefix, "_edges.tsv.gz")
  )
  expect_true(all(file.exists(expected_files)))

  message("Gene-gene graph test passed.")
  invisible(result)
}

main <- function() {
  test_that("gene-gene graph preparation works", {
    result <- test_prepare_gene_gene_graph()
    expect_s3_class(result$graph, "igraph")
  })
}

if (sys.nframe() == 0) {
  main()
}
