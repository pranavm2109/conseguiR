#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(testthat)
})

source("scripts/Internals/R/01_prepare_gene_reg_graph.R")

encode_reg_elements_path <- "data/raw/ENCODE/GRCh38-cCREs.bed"
encode_gene_links_zip_path <- "data/raw/ENCODE/Human-Gene-Links.zip"
default_gene_loc_path <- "data/raw/NCBI38/NCBI38.gene.loc"
default_gene_reg_test_output_dir <- "data/processed/test_outputs/gene_reg_graph"

make_gene_reg_test_path <- function(stem, ext = "") {
  dir.create(default_gene_reg_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_gene_reg_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

make_gene_reg_test_inputs <- function() {
  linked_ccres <- fread(
    cmd = sprintf(
      "unzip -p %s %s",
      shQuote(encode_gene_links_zip_path),
      shQuote(default_config$gene_link_members[[1]])
    ),
    header = FALSE,
    sep = "\t",
    fill = Inf,
    quote = "",
    nrows = 3000,
    showProgress = FALSE
  )
  reg_ids <- unique(as.character(linked_ccres$V1))
  reg_elements <- fread(encode_reg_elements_path, showProgress = FALSE)
  reg_elements <- reg_elements[as.character(V5) %in% reg_ids]
  reg_subset_path <- make_gene_reg_test_path("encode_ccres_subset", ".bed")
  fwrite(reg_elements, reg_subset_path, sep = "\t", col.names = FALSE)

  list(
    reg_elements_path = reg_subset_path,
    gene_links_zip_path = encode_gene_links_zip_path,
    gene_loc_path = default_gene_loc_path,
    reg_subset_ids = reg_ids
  )
}

subset_gene_links_zip <- function(zip_path, members, reg_ids, output_zip_path) {
  tmp_dir <- make_gene_reg_test_path("encode_gene_links_dir")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  output_zip_path_abs <- if (grepl("^/", output_zip_path)) {
    output_zip_path
  } else {
    file.path(repo_root, output_zip_path)
  }

  for (member in members) {
    dt <- fread(
      cmd = sprintf("unzip -p %s %s", shQuote(zip_path), shQuote(member)),
      header = FALSE,
      sep = "\t",
      fill = Inf,
      quote = "",
      showProgress = FALSE
    )
    dt <- dt[as.character(V1) %in% reg_ids]
    fwrite(dt, file.path(tmp_dir, member), sep = "\t", col.names = FALSE)
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  utils::zip(zipfile = output_zip_path_abs, files = members)
  normalizePath(output_zip_path_abs, winslash = "/", mustWork = TRUE)
}

test_prepare_gene_reg_graph <- function() {
  inputs <- make_gene_reg_test_inputs()
  output_prefix <- make_gene_reg_test_path("gene_reg_graph_test")
  gene_links_zip_subset <- make_gene_reg_test_path("gene_links_subset", ".zip")
  gene_links_zip_subset <- subset_gene_links_zip(
    zip_path = inputs$gene_links_zip_path,
    members = default_config$gene_link_members,
    reg_ids = inputs$reg_subset_ids,
    output_zip_path = gene_links_zip_subset
  )

  result <- prepare_gene_reg_graph(config = list(
    reg_elements_path = inputs$reg_elements_path,
    gene_links_zip_path = gene_links_zip_subset,
    gene_loc_path = inputs$gene_loc_path,
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
