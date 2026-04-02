#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

default_gene_reg_scoring_config <- list(
  graph_rds_path = "data/processed/gene_reg_graph_no_scores.rds",
  output_prefix = "data/processed/gene_reg_graph_scored"
)

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_gene_reg_graph_no_scores <- function(path = default_gene_reg_scoring_config$graph_rds_path) {
  if (!file.exists(path)) {
    stop("Gene-reg graph file does not exist: ", path)
  }

  graph <- readRDS(path)

  if (!inherits(graph, "igraph")) {
    stop("Expected an igraph object at: ", path)
  }

  graph
}

extract_gene_reg_graph_nodes <- function(graph) {
  if (!inherits(graph, "igraph")) {
    stop("`graph` must be an igraph object.")
  }

  as.data.table(as_data_frame(graph, what = "vertices"))
}

extract_gene_reg_graph_edges <- function(graph) {
  if (!inherits(graph, "igraph")) {
    stop("`graph` must be an igraph object.")
  }

  as.data.table(as_data_frame(graph, what = "edges"))
}

validate_gene_reg_graph_nodes <- function(nodes) {
  dt <- as.data.table(nodes)

  required_cols <- c("name", "node_id", "node_type")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Gene-reg node table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  bad_types <- setdiff(unique(dt$node_type), c("gene", "reg"))
  if (length(bad_types) > 0L) {
    stop("Gene-reg node table contains unsupported `node_type` values: ", paste(bad_types, collapse = ", "))
  }

  dt
}

standardize_gene_score_table <- function(scores, modality_name) {
  dt <- as.data.table(scores)

  if (!"gene_id" %in% names(dt) && "feature_id" %in% names(dt)) {
    dt[, gene_id := feature_id]
  }

  required_cols <- c("gene_id", "zstat")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      "Gene ",
      modality_name,
      " score table is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  unique(dt[, .(
    gene_id = toupper(as.character(gene_id)),
    zstat = as.numeric(zstat)
  )])[!is.na(gene_id) & gene_id != ""]
}

standardize_reg_score_table <- function(scores, modality_name) {
  dt <- as.data.table(scores)

  if (!"reg_elem_id" %in% names(dt) && "feature_id" %in% names(dt)) {
    dt[, reg_elem_id := feature_id]
  }

  required_cols <- c("reg_elem_id", "zstat")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      "Regulatory-element ",
      modality_name,
      " score table is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  unique(dt[, .(
    reg_elem_id = as.character(reg_elem_id),
    zstat = as.numeric(zstat)
  )])[!is.na(reg_elem_id) & reg_elem_id != ""]
}

merge_gene_modality_scores <- function(nodes, scores, score_col, modality_name) {
  dt <- copy(nodes)
  gene_nodes <- dt[node_type == "gene"]

  if (is.null(scores)) {
    gene_nodes[, (score_col) := 0]
    return(rbindlist(list(gene_nodes, dt[node_type != "gene"]), use.names = TRUE, fill = TRUE))
  }

  score_dt <- standardize_gene_score_table(scores, modality_name)
  gene_nodes[, gene_id := toupper(as.character(node_id))]
  gene_nodes <- merge(
    gene_nodes,
    score_dt,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
  gene_nodes[, (score_col) := fifelse(is.na(zstat), 0, zstat)]
  gene_nodes[, c("gene_id", "zstat") := NULL]

  rbindlist(list(gene_nodes, dt[node_type != "gene"]), use.names = TRUE, fill = TRUE)
}

merge_reg_modality_scores <- function(nodes, scores, score_col, modality_name) {
  dt <- copy(nodes)
  reg_nodes <- dt[node_type == "reg"]

  if (is.null(scores)) {
    reg_nodes[, (score_col) := 0]
    return(rbindlist(list(dt[node_type != "reg"], reg_nodes), use.names = TRUE, fill = TRUE))
  }

  score_dt <- standardize_reg_score_table(scores, modality_name)
  reg_nodes[, reg_elem_id := as.character(node_id)]
  reg_nodes <- merge(
    reg_nodes,
    score_dt,
    by = "reg_elem_id",
    all.x = TRUE,
    sort = FALSE
  )
  reg_nodes[, (score_col) := fifelse(is.na(zstat), 0, zstat)]
  reg_nodes[, c("reg_elem_id", "zstat") := NULL]

  rbindlist(list(dt[node_type != "reg"], reg_nodes), use.names = TRUE, fill = TRUE)
}

impose_scores_on_gene_reg_nodes <- function(
  nodes,
  gene_somatic_scores = NULL,
  gene_germline_scores = NULL,
  reg_somatic_scores = NULL,
  reg_germline_scores = NULL,
  reg_epigenomic_scores = NULL
) {
  dt <- validate_gene_reg_graph_nodes(nodes)

  dt <- merge_gene_modality_scores(dt, gene_somatic_scores, "somatic_score", "somatic")
  dt <- merge_gene_modality_scores(dt, gene_germline_scores, "germline_score", "germline")
  dt <- merge_reg_modality_scores(dt, reg_somatic_scores, "somatic_score", "somatic")
  dt <- merge_reg_modality_scores(dt, reg_germline_scores, "germline_score", "germline")
  dt <- merge_reg_modality_scores(dt, reg_epigenomic_scores, "epigenomic_score", "epigenomic")

  if (!"epigenomic_score" %in% names(dt)) {
    dt[, epigenomic_score := 0]
  }

  dt[node_type == "gene", epigenomic_score := 0]

  dt[, somatic_score := fifelse(is.na(somatic_score), 0, somatic_score)]
  dt[, germline_score := fifelse(is.na(germline_score), 0, germline_score)]
  dt[, epigenomic_score := fifelse(is.na(epigenomic_score), 0, epigenomic_score)]

  setcolorder(
    dt,
    c(
      "name",
      "node_id",
      "node_type",
      "somatic_score",
      "germline_score",
      "epigenomic_score",
      setdiff(names(dt), c("name", "node_id", "node_type", "somatic_score", "germline_score", "epigenomic_score"))
    )
  )

  setorder(dt, node_type, node_id)
  dt
}

impose_scores_on_gene_reg_graph <- function(
  graph,
  gene_somatic_scores = NULL,
  gene_germline_scores = NULL,
  reg_somatic_scores = NULL,
  reg_germline_scores = NULL,
  reg_epigenomic_scores = NULL
) {
  nodes <- extract_gene_reg_graph_nodes(graph)

  scored_nodes <- impose_scores_on_gene_reg_nodes(
    nodes = nodes,
    gene_somatic_scores = gene_somatic_scores,
    gene_germline_scores = gene_germline_scores,
    reg_somatic_scores = reg_somatic_scores,
    reg_germline_scores = reg_germline_scores,
    reg_epigenomic_scores = reg_epigenomic_scores
  )

  vertex_df <- as.data.frame(scored_nodes)
  edge_df <- as.data.frame(extract_gene_reg_graph_edges(graph))

  scored_graph <- graph_from_data_frame(
    d = edge_df,
    vertices = vertex_df,
    directed = is_directed(graph)
  )

  list(
    graph = scored_graph,
    nodes = scored_nodes,
    edges = as.data.table(edge_df)
  )
}

save_scored_gene_reg_graph_outputs <- function(graph, nodes, edges, output_prefix) {
  ensure_parent_dir(output_prefix)

  saveRDS(graph, paste0(output_prefix, ".rds"))
  fwrite(nodes, paste0(output_prefix, "_nodes.tsv.gz"), sep = "\t")
  fwrite(edges, paste0(output_prefix, "_edges.tsv.gz"), sep = "\t")

  invisible(output_prefix)
}

prepare_scored_gene_reg_graph <- function(
  graph = NULL,
  graph_rds_path = default_gene_reg_scoring_config$graph_rds_path,
  output_prefix = default_gene_reg_scoring_config$output_prefix,
  gene_somatic_scores = NULL,
  gene_germline_scores = NULL,
  reg_somatic_scores = NULL,
  reg_germline_scores = NULL,
  reg_epigenomic_scores = NULL,
  save_outputs = TRUE
) {
  if (is.null(graph)) {
    graph <- read_gene_reg_graph_no_scores(graph_rds_path)
  }

  result <- impose_scores_on_gene_reg_graph(
    graph = graph,
    gene_somatic_scores = gene_somatic_scores,
    gene_germline_scores = gene_germline_scores,
    reg_somatic_scores = reg_somatic_scores,
    reg_germline_scores = reg_germline_scores,
    reg_epigenomic_scores = reg_epigenomic_scores
  )

  if (isTRUE(save_outputs)) {
    save_scored_gene_reg_graph_outputs(
      graph = result$graph,
      nodes = result$nodes,
      edges = result$edges,
      output_prefix = output_prefix
    )
  }

  result
}
