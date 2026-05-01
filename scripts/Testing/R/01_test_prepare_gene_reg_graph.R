#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(testthat)
})

source("scripts/Internals/R/01_prepare_gene_reg_graph.R")

fixture_dir <- "scripts/Testing/fixtures/encode_gene_reg"
encode_reg_elements_path <- file.path(fixture_dir, "GRCh38-cCREs-tiny.bed")
default_gene_loc_path <- file.path(fixture_dir, "NCBI38.gene.loc.tiny")
default_gene_link_members <- c(
  "V4-hg38.Gene-Links.3D-Chromatin.txt",
  "V4-hg38.Gene-Links.CRISPR.txt",
  "V4-hg38.Gene-Links.eQTLs.txt"
)
default_gene_reg_test_output_dir <- file.path(tempdir(), "conseguiR_test_outputs", "gene_reg_graph")

make_gene_reg_test_path <- function(stem, ext = "") {
  dir.create(default_gene_reg_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_gene_reg_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

build_fixture_gene_links_zip <- function(members, output_zip_path) {
  tmp_dir <- make_gene_reg_test_path("encode_gene_links_dir")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  output_zip_path_abs <- if (grepl("^/", output_zip_path)) {
    output_zip_path
  } else {
    file.path(repo_root, output_zip_path)
  }

  for (member in members) {
    file.copy(
      from = file.path(fixture_dir, member),
      to = file.path(tmp_dir, member),
      overwrite = TRUE
    )
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(tmp_dir)
  zip_cmd <- getOption("zip")
  zip_available <- is.character(zip_cmd) &&
    length(zip_cmd) == 1L &&
    nzchar(zip_cmd) &&
    nzchar(Sys.which(zip_cmd))

  if (zip_available) {
    utils::zip(zipfile = output_zip_path_abs, files = members)
  } else {
    python_cmd <- Sys.which("python3")
    if (!nzchar(python_cmd)) {
      python_cmd <- Sys.which("python")
    }
    if (!nzchar(python_cmd)) {
      stop("Neither an external `zip` command nor Python is available to create fixture archives.")
    }

    zip_args <- c("-m", "zipfile", "-c", output_zip_path_abs, members)
    status <- system2(python_cmd, args = zip_args)
    if (!identical(status, 0L)) {
      stop("Failed to create fixture gene-link zip with Python fallback.")
    }
  }
  normalizePath(output_zip_path_abs, winslash = "/", mustWork = TRUE)
}

test_prepare_gene_reg_graph <- function() {
  output_prefix <- make_gene_reg_test_path("gene_reg_graph_test")
  gene_links_zip_subset <- make_gene_reg_test_path("gene_links_subset", ".zip")
  gene_links_zip_subset <- build_fixture_gene_links_zip(
    members = default_gene_link_members,
    output_zip_path = gene_links_zip_subset
  )

  result <- prepare_gene_reg_graph(config = list(
    reg_elements_path = encode_reg_elements_path,
    gene_links_zip_path = gene_links_zip_subset,
    gene_loc_path = default_gene_loc_path,
    gene_link_members = default_gene_link_members,
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
