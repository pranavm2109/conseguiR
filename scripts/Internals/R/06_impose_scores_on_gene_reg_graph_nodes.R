#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

default_gene_reg_scoring_config <- list(
  graph_rds_path = "data/processed/gene_reg_graph_no_scores.rds",
  output_prefix = "data/processed/gene_reg_graph_scored",
  filter_to_supported_universe = FALSE
)

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_gene_reg_graph_no_scores <- function(path = default_gene_reg_scoring_config$graph_rds_path) {
  if (file.exists(path)) {
    graph <- readRDS(path)

    if (!inherits(graph, "igraph")) {
      stop("Expected an igraph object at: ", path)
    }

    return(graph)
  }

  nodes_path <- sub("\\.rds$", "_nodes.tsv.gz", path)
  edges_path <- sub("\\.rds$", "_edges.tsv.gz", path)
  if (!identical(nodes_path, path) && file.exists(nodes_path) && file.exists(edges_path)) {
    nodes <- data.table::as.data.table(data.table::fread(nodes_path, showProgress = FALSE))
    edges <- data.table::as.data.table(data.table::fread(edges_path, showProgress = FALSE))

    if (!all(c("node_id", "node_type") %in% names(nodes))) {
      stop("Gene-reg node table is missing required columns needed to reconstruct the graph.")
    }
    if (!all(c("from", "to") %in% names(edges))) {
      stop("Gene-reg edge table is missing required columns needed to reconstruct the graph.")
    }

    vertices <- as.data.frame(nodes)
    if ("name" %in% names(vertices)) {
      vertices$name <- as.character(vertices$name)
    } else {
      vertices$name <- as.character(vertices$node_id)
    }

    return(graph_from_data_frame(
      d = as.data.frame(edges),
      vertices = vertices,
      directed = TRUE
    ))
  }

  stop(
    "Gene-reg graph file does not exist: ", path,
    ". Neither a seed `.rds` graph nor matching `_nodes.tsv.gz`/`_edges.tsv.gz` files were found."
  )
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

read_backend_gene_id_map <- function() {
  helper <- tryCatch(
    getFromNamespace(".conseguiR_backend_resource_path", "conseguiR"),
    error = function(e) NULL
  )
  loc_path <- if (!is.null(helper)) {
    helper("NCBI38.gene.loc")
  } else {
    loc_candidates <- c(
      getOption("conseguiR.backend_resource_dir", NULL),
      file.path(tools::R_user_dir("conseguiR", which = "cache"), "backend_resources"),
      "inst/extdata/backend",
      "extdata/backend",
      "data/raw/NCBI38/NCBI38.gene.loc"
    )
    loc_candidates <- unique(as.character(loc_candidates))
    loc_candidates <- loc_candidates[!is.na(loc_candidates) & nzchar(loc_candidates)]
    if (!any(grepl("NCBI38\\.gene\\.loc$", loc_candidates))) {
      loc_candidates <- file.path(loc_candidates, "NCBI38.gene.loc")
    }
    existing <- loc_candidates[file.exists(loc_candidates)]
    if (length(existing) == 0L) NULL else existing[[1]]
  }

  if (is.null(loc_path) || !nzchar(loc_path)) {
    stop("Could not locate backend gene ID map file: NCBI38.gene.loc")
  }

  dt <- data.table::fread(loc_path, header = FALSE, showProgress = FALSE)
  if (ncol(dt) < 6L) {
    stop("Backend gene ID map does not contain the expected symbol column.")
  }

  unique(dt[, .(
    feature_id = as.character(V1),
    gene_id = toupper(as.character(V6))
  )], by = "feature_id")[!is.na(feature_id) & feature_id != "" & !is.na(gene_id) & gene_id != ""]
}

standardize_gene_score_table <- function(scores, modality_name) {
  dt <- as.data.table(scores)

  if (!"gene_id" %in% names(dt) && "feature_id" %in% names(dt)) {
    id_map <- read_backend_gene_id_map()
    dt[, feature_id := as.character(feature_id)]
    dt <- merge(dt, id_map, by = "feature_id", all.x = TRUE, sort = FALSE)
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

  dt <- dt[, .(
    gene_id = toupper(as.character(gene_id)),
    zstat = as.numeric(zstat)
  )][!is.na(gene_id) & gene_id != "" & !is.na(zstat)]

  # If duplicate rows are supplied for the same feature, keep the most
  # vulnerability-aligned entry rather than the largest-magnitude signed score.
  dt <- dt[order(-zstat)][, .SD[1], by = gene_id]

  dt
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

  dt <- dt[, .(
    reg_elem_id = as.character(reg_elem_id),
    zstat = as.numeric(zstat)
  )][!is.na(reg_elem_id) & reg_elem_id != "" & !is.na(zstat)]

  dt <- dt[order(-zstat)][, .SD[1], by = reg_elem_id]

  dt
}

merge_gene_modality_scores <- function(nodes, scores, score_col, modality_name) {
  dt <- data.table::copy(nodes)
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
  dt <- data.table::copy(nodes)
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

compute_supported_gene_ids <- function(gene_somatic_scores = NULL, gene_germline_scores = NULL) {
  gene_sets <- list()

  if (!is.null(gene_somatic_scores)) {
    gene_sets <- c(gene_sets, list(unique(standardize_gene_score_table(gene_somatic_scores, "somatic")$gene_id)))
  }
  if (!is.null(gene_germline_scores)) {
    gene_sets <- c(gene_sets, list(unique(standardize_gene_score_table(gene_germline_scores, "germline")$gene_id)))
  }

  if (length(gene_sets) == 0L) {
    return(NULL)
  }

  Reduce(intersect, gene_sets)
}

compute_supported_reg_ids <- function(
  reg_somatic_scores = NULL,
  reg_germline_scores = NULL,
  reg_epigenomic_scores = NULL
) {
  if (!is.null(reg_epigenomic_scores)) {
    return(unique(standardize_reg_score_table(reg_epigenomic_scores, "epigenomic")$reg_elem_id))
  }

  reg_sets <- list()
  if (!is.null(reg_somatic_scores)) {
    reg_sets <- c(reg_sets, list(unique(standardize_reg_score_table(reg_somatic_scores, "somatic")$reg_elem_id)))
  }
  if (!is.null(reg_germline_scores)) {
    reg_sets <- c(reg_sets, list(unique(standardize_reg_score_table(reg_germline_scores, "germline")$reg_elem_id)))
  }

  if (length(reg_sets) == 0L) {
    return(NULL)
  }

  Reduce(intersect, reg_sets)
}

subset_gene_reg_graph_to_supported_universe <- function(
  graph,
  gene_somatic_scores = NULL,
  gene_germline_scores = NULL,
  reg_somatic_scores = NULL,
  reg_germline_scores = NULL,
  reg_epigenomic_scores = NULL
) {
  nodes <- validate_gene_reg_graph_nodes(extract_gene_reg_graph_nodes(graph))
  edges <- extract_gene_reg_graph_edges(graph)

  keep_gene_ids <- compute_supported_gene_ids(
    gene_somatic_scores = gene_somatic_scores,
    gene_germline_scores = gene_germline_scores
  )
  keep_reg_ids <- compute_supported_reg_ids(
    reg_somatic_scores = reg_somatic_scores,
    reg_germline_scores = reg_germline_scores,
    reg_epigenomic_scores = reg_epigenomic_scores
  )

  if (is.null(keep_gene_ids) || length(keep_gene_ids) == 0L) {
    stop("Could not derive a supported gene universe from the supplied pre-diffusion gene score tables.")
  }
  if (is.null(keep_reg_ids) || length(keep_reg_ids) == 0L) {
    stop("Could not derive a supported regulatory universe from the supplied regulatory score tables.")
  }

  edge_dt <- as.data.table(edges)
  edge_dt <- edge_dt[toupper(from) %in% keep_gene_ids & to %in% keep_reg_ids]

  if (nrow(edge_dt) == 0L) {
    stop("Supported-universe filtering removed all gene-reg edges. Check identifier harmonization and regulatory coverage.")
  }

  kept_gene_ids <- unique(toupper(edge_dt$from))
  kept_reg_ids <- unique(edge_dt$to)

  node_dt <- as.data.table(nodes)[
    (node_type == "gene" & toupper(node_id) %in% kept_gene_ids) |
      (node_type == "reg" & node_id %in% kept_reg_ids)
  ]

  if (nrow(node_dt) == 0L) {
    stop("Supported-universe filtering removed all graph nodes.")
  }

  graph_from_data_frame(
    d = as.data.frame(edge_dt),
    vertices = as.data.frame(node_dt),
    directed = is_directed(graph)
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
  filter_to_supported_universe = default_gene_reg_scoring_config$filter_to_supported_universe,
  save_outputs = TRUE
) {
  if (is.null(graph)) {
    graph <- read_gene_reg_graph_no_scores(graph_rds_path)
  }

  if (isTRUE(filter_to_supported_universe)) {
    graph <- subset_gene_reg_graph_to_supported_universe(
      graph = graph,
      gene_somatic_scores = gene_somatic_scores,
      gene_germline_scores = gene_germline_scores,
      reg_somatic_scores = reg_somatic_scores,
      reg_germline_scores = reg_germline_scores,
      reg_epigenomic_scores = reg_epigenomic_scores
    )
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
