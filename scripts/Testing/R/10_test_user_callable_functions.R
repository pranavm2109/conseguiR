#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
  library(devtools)
})

source("R/backend_resources.R")
source("R/user_api.R")
source("scripts/Internals/R/06_impose_scores_on_gene_reg_graph_nodes.R")
source("R/zzz.R")
# Manually trigger package initialization before sourcing the external wrapper layer.
.onLoad(libname = ".", pkgname = "conseguiR")
devtools::load_all(".", quiet = TRUE)

source("scripts/Externals/R/00_user_callable_functions.R")

default_external_test_output_dir <- file.path(tempdir(), "conseguiR_test_outputs", "external_api")
default_external_gwas_path <- "data/raw/Testing/34737426-GCST90043906-EFO_0000403.h.tsv"
default_external_somatic_path <- "data/raw/Testing/2026-01-09_no_CLL_lymph_only_pcawg_maf_tcga_order_hg38.maf"
default_external_reg_ref_path <- "data/raw/Testing/reg_elements_valid.loc"
default_external_bw_files <- c(
  "data/raw/Testing/SRR1020514_DLBCL_P265_H3K27ac_ChIPseq.bw",
  "data/raw/Testing/SRR1020516_DLBCL_P286_H3K27ac_ChIPseq.bw",
  "data/raw/Testing/SRR1020518_DLBCL_P397_H3K27ac_ChIPseq.bw"
)

make_live_validation_inputs <- function() {
  list(
    gwas = fread(default_external_gwas_path, nrows = 2000L, showProgress = FALSE),
    somatic = fread(default_external_somatic_path, nrows = 2000L, showProgress = FALSE)
  )
}

make_external_test_path <- function(stem, ext = "") {
  dir.create(default_external_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_external_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

test_backend_dir <- file.path(tempdir(), "conseguiR_backend_test_10")
options(conseguiR.backend_dir = test_backend_dir)
.conseguiR_initialize_backend_graphs(
  backend_dir = test_backend_dir,
  build_gene_reg = TRUE,
  build_gene_gene = TRUE,
  force = TRUE,
  strict = FALSE,
  quiet = TRUE
)

make_existing_score_bundles <- function() {
  graph <- read_gene_reg_graph_no_scores(.conseguiR_backend_paths(test_backend_dir)$gene_reg_graph_rds)
  nodes <- extract_gene_reg_graph_nodes(graph)
  gene_ids <- unique(nodes[node_type == "gene", node_id])[1:6]
  reg_ids <- unique(nodes[node_type == "reg", node_id])[1:6]

  germline <- new_bundle(
    type = "germline_scores",
    objects = list(
      gene_scores = data.table(gene_id = gene_ids, zstat = c(2.2, 1.7, 0.9, -0.4, 0.3, -0.2)),
      reg_scores = data.table(reg_elem_id = reg_ids, zstat = c(1.8, 0.7, -0.6, 0.4, 0.2, -0.1))
    )
  )

  somatic <- new_bundle(
    type = "somatic_scores",
    objects = list(
      gene_scores = data.table(gene_id = gene_ids, zstat = c(1.9, -1.1, 0.8, 0.5, -0.2, 0.1)),
      reg_scores = data.table(reg_elem_id = reg_ids, zstat = c(1.4, -0.5, 0.9, 0.3, -0.1, 0.2))
    )
  )

  epigenomic <- new_bundle(
    type = "epigenomic_scores",
    objects = list(
      reg_scores = data.table(reg_elem_id = reg_ids, zstat = c(2.8, -1.5, 1.0, 0.6, -0.3, 0.4))
    )
  )

  list(
    germline = germline,
    somatic = somatic,
    epigenomic = epigenomic
  )
}

materialize_plaintext_gene_gene_paths <- function() {
  backend_paths <- .conseguiR_backend_paths(test_backend_dir)
  if (!file.exists(backend_paths$gene_gene_graph_rds)) {
    skip("Gene-gene backend graph is unavailable in this environment.")
  }
  graph <- readRDS(backend_paths$gene_gene_graph_rds)

  nodes_path <- make_external_test_path("gene_gene_nodes_plain", ".tsv")
  edges_path <- make_external_test_path("gene_gene_edges_plain", ".tsv")

  fwrite(as.data.table(igraph::as_data_frame(graph, what = "vertices")), nodes_path, sep = "\t")
  fwrite(as.data.table(igraph::as_data_frame(graph, what = "edges")), edges_path, sep = "\t")

  list(
    gg_nodes_path = nodes_path,
    gg_edges_path = edges_path
  )
}

test_validate_inputs_external_live <- function() {
  inputs <- make_live_validation_inputs()

  result <- validate_inputs(
    gwas_sumstats = inputs$gwas,
    somatic_maf = inputs$somatic,
    reg_ref_path = default_external_reg_ref_path,
    epigenomic_tracks = default_external_bw_files
  )

  expect_s3_class(result, "conseguiR_bundle")
  expect_identical(result$bundle_type, "validation")
  expect_true(is.data.frame(result$gwas))
  expect_true(is.data.frame(result$somatic_maf))
  expect_true(inherits(result$regulatory_elements, "GRanges"))
  expect_true(is.list(result$epigenomic))
  expect_true(isTRUE(result$config$has_gwas))
  expect_true(isTRUE(result$config$has_somatic_maf))
  expect_true(isTRUE(result$config$has_reg_ref))
  expect_true(isTRUE(result$config$has_epigenomic))
}

test_external_graph_to_plot_chain_live <- function(print_outputs = TRUE) {
  bundles <- make_existing_score_bundles()

  scored_prefix <- make_external_test_path("external_scored_graph")
  scored <- build_scored_gene_reg_graph(
    output_prefix = scored_prefix,
    germline_scores = bundles$germline,
    somatic_scores = bundles$somatic,
    epigenomic_scores = bundles$epigenomic,
    save_outputs = TRUE
  )

  expect_s3_class(scored, "conseguiR_bundle")
  expect_true(file.exists(paste0(scored_prefix, ".rds")))
  expect_true(file.exists(paste0(scored_prefix, "_nodes.tsv.gz")))
  expect_true(file.exists(paste0(scored_prefix, "_edges.tsv.gz")))
  expect_true(is.data.frame(scored$nodes))
  expect_true(is.data.frame(scored$edges))

  diffusion_dir <- make_external_test_path("external_diffusion_dir")
  dir.create(diffusion_dir, recursive = TRUE, showWarnings = FALSE)
  diffusion <- run_gene_reg_diffusion(
    scored_graph = scored,
    output_dir = diffusion_dir,
    output_stem = "external_diffusion",
    top_n_to_save = 15L
  )

  expect_s3_class(diffusion, "conseguiR_bundle")
  expect_true(file.exists(diffusion$output_paths$all_genes_path))
  expect_true(file.exists(diffusion$output_paths$top_genes_path))
  expect_true(is.data.frame(diffusion$all_genes))
  expect_true(is.data.frame(diffusion$top_genes))

  subgraph_dir <- make_external_test_path("external_subgraph_dir")
  dir.create(subgraph_dir, recursive = TRUE, showWarnings = FALSE)
  gene_gene_paths <- materialize_plaintext_gene_gene_paths()
  selected <- call_selected_subgraph(
    diffusion = diffusion,
    gg_nodes_path = gene_gene_paths$gg_nodes_path,
    gg_edges_path = gene_gene_paths$gg_edges_path,
    output_dir = subgraph_dir,
    output_stem = "external_selected_subgraph",
    target_genes = 12L,
    candidate_pool_size = 60L,
    max_time_seconds = 20L,
    num_workers = 2L
  )

  expect_s3_class(selected, "conseguiR_bundle")
  expect_true(file.exists(selected$output_paths$nodes_path))
  expect_true(file.exists(selected$output_paths$edges_path))
  expect_true(file.exists(selected$output_paths$summary_path))
  expect_true(file.exists(selected$output_paths$graphml_path))
  expect_true(is.data.frame(selected$nodes))
  expect_true(is.data.frame(selected$edges))
  expect_true(is.data.frame(selected$summary))

  plot_prefix <- make_external_test_path("external_plot_bundle")
  plot_path <- make_external_test_path("external_selected_subgraph", ".pdf")
  plotted <- plot_selected_subgraph(
    selected_subgraph = selected,
    bundle_output_prefix = plot_prefix,
    plot_file_path = plot_path,
    title = "External API Test Plot",
    top_n_labels = Inf,
    save_bundle = TRUE,
    save_plot = TRUE
  )

  expect_s3_class(plotted, "conseguiR_bundle")
  expect_true(file.exists(plot_path))
  expect_true(file.exists(paste0(plot_prefix, ".rds")))
  expect_true(file.exists(paste0(plot_prefix, "_nodes.tsv.gz")))
  expect_true(file.exists(paste0(plot_prefix, "_edges.tsv.gz")))
  expect_true(inherits(plotted$plot, "ggplot"))

  if (isTRUE(print_outputs)) {
    message("Top externally selected genes:")
    print(selected$nodes[order(-prize)][1:min(10L, .N), .(gene_name, prize)])
    message("Saved external API test plot to: ", plot_path)
  }

  invisible(list(
    scored = scored,
    diffusion = diffusion,
    selected = selected,
    plotted = plotted
  ))
}

test_validate_inputs_external_negative_missing_reg_ref_for_epigenomic <- function() {
  expect_error(
    validate_inputs(
      epigenomic_tracks = default_external_bw_files
    ),
    regexp = "`reg_ref_path` is required when validating epigenomic tracks.",
    fixed = TRUE
  )
}

test_build_scored_gene_reg_graph_external_negative_bad_scores <- function() {
  expect_error(
    build_scored_gene_reg_graph(
      gene_germline_scores = data.table(gene = "MYC", value = 1.2),
      save_outputs = FALSE
    ),
    regexp = "Gene germline score table is missing required columns",
    fixed = TRUE
  )
}

test_run_gene_reg_diffusion_external_negative_missing_paths <- function() {
  expect_error(
    run_gene_reg_diffusion(
      nodes_path = "data/processed/does_not_exist_nodes.tsv.gz",
      edges_path = "data/processed/does_not_exist_edges.tsv.gz",
      output_dir = make_external_test_path("bad_diffusion_dir"),
      output_stem = "bad_diffusion"
    ),
    regexp = "Scored gene-reg node file does not exist",
    fixed = TRUE
  )
}

test_call_selected_subgraph_external_negative_missing_diffusion <- function() {
  expect_error(
    call_selected_subgraph(
      diffusion_path = "data/processed/does_not_exist_diffusion.tsv",
      output_dir = make_external_test_path("bad_subgraph_dir"),
      output_stem = "bad_subgraph"
    ),
    regexp = "Diffusion results file does not exist",
    fixed = TRUE
  )
}

test_plot_selected_subgraph_external_negative_missing_plot_path <- function() {
  fixture <- new_bundle(
    type = "selected_subgraph",
    objects = list(
      nodes = data.table(
        node_id = c("G1", "G2"),
        gene_name = c("GENE1", "GENE2"),
        prize = c(2, 1),
        post_integrated = c(2.0, 1.0),
        post_vulnerability = c(2.1, 1.1),
        post_norm = c(2.2, 1.2)
      ),
      edges = data.table(
        gene_u = "G1",
        gene_v = "G2"
      ),
      summary = data.table(
        solver_status = "OPTIMAL",
        target_genes = 2L,
        n_selected_nodes = 2L,
        n_selected_edges = 1L
      )
    )
  )

  expect_error(
    plot_selected_subgraph(
      selected_subgraph = fixture,
      save_plot = TRUE,
      plot_file_path = ""
    ),
    regexp = "`plot_file_path` must be provided when `save_plot = TRUE`.",
    fixed = TRUE
  )
}

main <- function() {
  test_that("external validate_inputs works on the live in-repo inputs", {
    test_validate_inputs_external_live()
  })

  test_that("external graph-to-plot chain runs on the live in-repo outputs", {
    test_external_graph_to_plot_chain_live(print_outputs = TRUE)
  })

  test_that("external validate_inputs fails clearly for epigenomic validation without a regulatory reference", {
    test_validate_inputs_external_negative_missing_reg_ref_for_epigenomic()
  })

  test_that("external scored graph building fails clearly for malformed score tables", {
    test_build_scored_gene_reg_graph_external_negative_bad_scores()
  })

  test_that("external diffusion call fails clearly when scored graph paths are missing", {
    test_run_gene_reg_diffusion_external_negative_missing_paths()
  })

  test_that("external subgraph calling fails clearly when diffusion input is missing", {
    test_call_selected_subgraph_external_negative_missing_diffusion()
  })

  test_that("external plotting fails clearly when a plot save path is missing", {
    test_plot_selected_subgraph_external_negative_missing_plot_path()
  })
}

if (sys.nframe() == 0) {
  main()
}
