#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

conseguiR_runtime_file <- function(relpath) {
  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(candidate)
  }

  pkg_path <- system.file(relpath, package = "conseguiR")
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }

  stop("Could not locate required runtime file: ", relpath)
}

default_subgraph_config <- list(
  diffusion_path = "data/processed/gene_reg_graph_diffusion_all_genes.tsv",
  gg_nodes_path = "data/processed/gene_gene_graph_nodes.tsv.gz",
  gg_edges_path = "data/processed/gene_gene_graph_edges.tsv.gz",
  output_dir = "data/processed",
  output_stem = "gene_gene_selected_subgraph",
  target_genes = 50L,
  candidate_pool_size = 400L,
  min_confidence = 0.0,
  max_edges_in_model = 12000L,
  node_prize_weight = 1.0,
  edge_conf_weight = 1.0,
  edge_cost_weight = 1.0,
  node_scale = 1000L,
  edge_scale = 1000L,
  max_time_seconds = 600L,
  num_workers = 8L,
  random_seed = 42L,
  prize_column = "post_norm",
  confidence_column = "confidence",
  edge_cost_column = "weight"
)

## Python execution is now delegated to the shared basilisk helper sourced by
## the user-facing external API runtime.

read_diffusion_results <- function(path = default_subgraph_config$diffusion_path) {
  if (!file.exists(path)) {
    stop("Diffusion results file does not exist: ", path)
  }

  data.table::as.data.table(data.table::fread(path, sep = "\t", showProgress = FALSE))
}

read_gene_gene_nodes <- function(path = default_subgraph_config$gg_nodes_path) {
  if (!file.exists(path)) {
    stop("Gene-gene node file does not exist: ", path)
  }

  as.data.table(fread(path, sep = "\t", showProgress = FALSE))
}

read_gene_gene_edges <- function(path = default_subgraph_config$gg_edges_path) {
  if (!file.exists(path)) {
    stop("Gene-gene edge file does not exist: ", path)
  }

  as.data.table(fread(path, sep = "\t", showProgress = FALSE))
}

validate_diffusion_results <- function(diffusion, prize_column = default_subgraph_config$prize_column) {
  dt <- as.data.table(diffusion)
  required_cols <- c("node_id", "gene_name")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Diffusion results are missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!(prize_column %in% names(dt))) {
    required_post_cols <- c("post_germline", "post_somatic", "post_epigenomic")
    missing_post_cols <- setdiff(required_post_cols, names(dt))
    if (length(missing_post_cols) > 0L) {
      stop(
        "Diffusion results are missing prize column `", prize_column,
        "` and do not contain all three post-diffusion component columns."
      )
    }
  }

  if (any(is.na(dt$node_id))) {
    stop("Diffusion results contain missing `node_id` values.")
  }

  if (any(is.na(dt$gene_name))) {
    stop("Diffusion results contain missing `gene_name` values.")
  }

  dt
}

validate_gene_gene_nodes <- function(nodes) {
  dt <- as.data.table(nodes)
  required_cols <- c("node_id", "node_type")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Gene-gene node table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(is.na(dt$node_id)) || any(duplicated(dt$node_id))) {
    stop("Gene-gene node table must have unique, non-missing node_id values.")
  }

  bad_types <- setdiff(unique(dt$node_type), "gene")
  if (length(bad_types) > 0L) {
    stop("Gene-gene node table contains unsupported node_type values: ", paste(bad_types, collapse = ", "))
  }

  dt
}

validate_gene_gene_edges <- function(edges, confidence_column = default_subgraph_config$confidence_column, edge_cost_column = default_subgraph_config$edge_cost_column) {
  dt <- as.data.table(edges)
  required_cols <- c("from", "to", confidence_column, edge_cost_column)
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Gene-gene edge table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(is.na(dt$from)) || any(is.na(dt$to))) {
    stop("Gene-gene edge table contains missing endpoints.")
  }

  if (any(is.na(dt[[confidence_column]]))) {
    stop("Gene-gene edge table contains missing `", confidence_column, "` values.")
  }

  if (any(is.na(dt[[edge_cost_column]]))) {
    stop("Gene-gene edge table contains missing `", edge_cost_column, "` values.")
  }

  dt
}

run_cardinality_subgraph_calling <- function(
  diffusion_path = default_subgraph_config$diffusion_path,
  gg_nodes_path = default_subgraph_config$gg_nodes_path,
  gg_edges_path = default_subgraph_config$gg_edges_path,
  output_dir = default_subgraph_config$output_dir,
  output_stem = default_subgraph_config$output_stem,
  target_genes = default_subgraph_config$target_genes,
  candidate_pool_size = default_subgraph_config$candidate_pool_size,
  min_confidence = default_subgraph_config$min_confidence,
  max_edges_in_model = default_subgraph_config$max_edges_in_model,
  node_prize_weight = default_subgraph_config$node_prize_weight,
  edge_conf_weight = default_subgraph_config$edge_conf_weight,
  edge_cost_weight = default_subgraph_config$edge_cost_weight,
  node_scale = default_subgraph_config$node_scale,
  edge_scale = default_subgraph_config$edge_scale,
  max_time_seconds = default_subgraph_config$max_time_seconds,
  num_workers = default_subgraph_config$num_workers,
  random_seed = default_subgraph_config$random_seed,
  prize_column = default_subgraph_config$prize_column,
  confidence_column = default_subgraph_config$confidence_column,
  edge_cost_column = default_subgraph_config$edge_cost_column,
  python_path = NULL
) {
  if (!is.null(python_path) && nzchar(python_path)) {
    warning(
      "`python_path` is deprecated and ignored. ",
      "conseguiR now runs selected-subgraph calling inside a managed basilisk environment.",
      call. = FALSE
    )
  }

  config <- list(
    diffusion_path = diffusion_path,
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path,
    output_dir = output_dir,
    output_stem = output_stem,
    target_genes = target_genes,
    candidate_pool_size = candidate_pool_size,
    min_confidence = min_confidence,
    max_edges_in_model = max_edges_in_model,
    node_prize_weight = node_prize_weight,
    edge_conf_weight = edge_conf_weight,
    edge_cost_weight = edge_cost_weight,
    node_scale = node_scale,
    edge_scale = edge_scale,
    max_time_seconds = max_time_seconds,
    num_workers = num_workers,
    random_seed = random_seed,
    prize_column = prize_column,
    confidence_column = confidence_column,
    edge_cost_column = edge_cost_column
  )

  output_paths <- .conseguiR_basilisk_run_python_module(
    script_path = conseguiR_runtime_file("scripts/Internals/Python/08_call_subgraph.py"),
    module_name = "conseguiR_step8_r",
    function_name = "run_cardinality_subgraph_calling",
    config = config,
    config_class_name = "SubgraphConfig"
  )

  list(
    config = config,
    output_paths = output_paths
  )
}
