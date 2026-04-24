#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

default_selected_subgraph_plot_config <- list(
  nodes_path = "data/processed/gene_gene_selected_subgraph_nodes.tsv",
  edges_path = "data/processed/gene_gene_selected_subgraph_edges.tsv",
  summary_path = "data/processed/gene_gene_selected_subgraph_summary.tsv",
  output_prefix = "data/processed/gene_gene_selected_subgraph_plot_bundle"
)

read_selected_subgraph_nodes <- function(path = default_selected_subgraph_plot_config$nodes_path) {
  if (!file.exists(path)) {
    stop("Selected subgraph node file does not exist: ", path)
  }

  as.data.table(fread(path))
}

read_selected_subgraph_edges <- function(path = default_selected_subgraph_plot_config$edges_path) {
  if (!file.exists(path)) {
    stop("Selected subgraph edge file does not exist: ", path)
  }

  as.data.table(fread(path))
}

read_selected_subgraph_summary <- function(path = default_selected_subgraph_plot_config$summary_path) {
  if (!file.exists(path)) {
    stop("Selected subgraph summary file does not exist: ", path)
  }

  as.data.table(fread(path))
}

validate_selected_subgraph_nodes <- function(nodes) {
  dt <- as.data.table(nodes)

  required_cols <- c("node_id", "gene_name", "prize")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Selected subgraph node table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(is.na(dt$node_id)) || any(dt$node_id == "")) {
    stop("Selected subgraph node table contains missing `node_id` values.")
  }

  if (anyDuplicated(dt$node_id)) {
    stop("Selected subgraph node table contains duplicated `node_id` values.")
  }

  if (any(is.na(dt$gene_name)) || any(dt$gene_name == "")) {
    stop("Selected subgraph node table contains missing `gene_name` values.")
  }

  if (!is.numeric(dt$prize)) {
    suppressWarnings(dt[, prize := as.numeric(prize)])
  }
  if (any(is.na(dt$prize))) {
    stop("Selected subgraph node table contains non-numeric `prize` values.")
  }

  dt
}

validate_selected_subgraph_edges <- function(edges) {
  dt <- as.data.table(edges)

  required_cols <- c("gene_u", "gene_v")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Selected subgraph edge table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(is.na(dt$gene_u)) || any(dt$gene_u == "") || any(is.na(dt$gene_v)) || any(dt$gene_v == "")) {
    stop("Selected subgraph edge table contains missing gene endpoints.")
  }

  dt
}

validate_selected_subgraph_summary <- function(summary) {
  dt <- as.data.table(summary)

  required_cols <- c("solver_status", "target_genes", "n_selected_nodes", "n_selected_edges")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Selected subgraph summary table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  dt
}

sanitize_nonfinite_numeric_columns <- function(dt) {
  out <- data.table::copy(as.data.table(dt))
  numeric_cols <- names(out)[vapply(out, is.numeric, logical(1))]

  for (col in numeric_cols) {
    values <- out[[col]]
    finite_mask <- is.finite(values)

    if (!all(finite_mask)) {
      if (any(finite_mask)) {
        finite_max <- max(values[finite_mask], na.rm = TRUE)
        values[!finite_mask] <- finite_max
      } else {
        values[!finite_mask] <- 0
      }
      out[[col]] <- values
    }
  }

  out
}

prepare_subgraph_plot_nodes <- function(nodes, sanitize_nonfinite = TRUE) {
  dt <- validate_selected_subgraph_nodes(nodes)

  if (isTRUE(sanitize_nonfinite)) {
    dt <- sanitize_nonfinite_numeric_columns(dt)
  }

  dt <- data.table::copy(dt)
  dt[, label := gene_name]
  dt[, node_class := "selected_gene"]
  dt
}

prepare_subgraph_plot_edges <- function(edges, nodes) {
  edge_dt <- validate_selected_subgraph_edges(edges)
  node_dt <- validate_selected_subgraph_nodes(nodes)

  valid_ids <- unique(node_dt$node_id)
  edge_dt <- edge_dt[gene_u %in% valid_ids & gene_v %in% valid_ids]

  if (nrow(edge_dt) == 0L) {
    stop("Selected subgraph edges do not overlap the selected node identifiers.")
  }

  edge_dt[, from := gene_u]
  edge_dt[, to := gene_v]
  edge_dt
}

build_selected_subgraph_igraph <- function(nodes, edges) {
  node_dt <- prepare_subgraph_plot_nodes(nodes)
  edge_dt <- prepare_subgraph_plot_edges(edges, node_dt)

  graph_from_data_frame(
    d = as.data.frame(edge_dt),
    vertices = as.data.frame(node_dt),
    directed = FALSE
  )
}

build_selected_subgraph_tidygraph <- function(nodes, edges) {
  if (!requireNamespace("tidygraph", quietly = TRUE)) {
    return(NULL)
  }

  tidygraph::as_tbl_graph(build_selected_subgraph_igraph(nodes, edges))
}

prepare_cytoscape_data_frames <- function(nodes, edges) {
  node_dt <- prepare_subgraph_plot_nodes(nodes)
  edge_dt <- prepare_subgraph_plot_edges(edges, node_dt)

  cy_nodes <- data.table::copy(node_dt)
  if (!"id" %in% names(cy_nodes)) {
    cy_nodes[, id := node_id]
  }
  cy_nodes[, shared.name := gene_name]
  cy_nodes[, name := gene_name]

  cy_edges <- data.table::copy(edge_dt)
  cy_edges[, source := from]
  cy_edges[, target := to]
  cy_edges[, interaction := "selected_subgraph_edge"]

  list(
    nodes = cy_nodes,
    edges = cy_edges
  )
}

create_selected_subgraph_cytoscape_network <- function(
  nodes,
  edges,
  title = "conseguiR Selected Gene Subgraph",
  collection = "conseguiR"
) {
  if (!requireNamespace("RCy3", quietly = TRUE)) {
    stop("RCy3 is required to create a Cytoscape network.")
  }

  if (!RCy3::cytoscapePing()) {
    stop("Cytoscape is not reachable from RCy3.")
  }

  cy_data <- prepare_cytoscape_data_frames(nodes, edges)

  RCy3::createNetworkFromDataFrames(
    nodes = as.data.frame(cy_data$nodes),
    edges = as.data.frame(cy_data$edges),
    title = title,
    collection = collection
  )
}

create_selected_subgraph_plot <- function(
  bundle,
  layout = "fr",
  top_n_labels = Inf,
  title = "conseguiR Selected Gene Subgraph",
  subtitle = NULL
) {
  required_packages <- c("ggplot2", "scales", "ggrepel")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, logical(1))]
  if (length(missing_packages) > 0L) {
    stop("The following packages are required to create the selected subgraph plot: ",
         paste(missing_packages, collapse = ", "))
  }

  graph_obj <- bundle$igraph
  if (is.null(graph_obj) || !inherits(graph_obj, "igraph")) {
    graph_obj <- build_selected_subgraph_igraph(bundle$nodes, bundle$edges)
  }

  node_dt <- as.data.table(bundle$nodes)
  if (is.infinite(top_n_labels)) {
    label_count <- nrow(node_dt)
  } else {
    label_count <- min(as.integer(top_n_labels), nrow(node_dt))
  }
  label_ids <- character()
  if (label_count > 0L) {
    label_ids <- node_dt[order(-prize, gene_name)][1:label_count, node_id]
  }

  label_text_size <- if (label_count > 35L) 2.8 else 3.4
  label_box_padding <- if (label_count > 35L) 0.22 else 0.35
  label_point_padding <- if (label_count > 35L) 0.08 else 0.15
  label_force <- if (label_count > 35L) 1.8 else 1.2

  coords <- switch(
    layout,
    fr = igraph::layout_with_fr(graph_obj),
    kk = igraph::layout_with_kk(graph_obj),
    lgl = igraph::layout_with_lgl(graph_obj),
    circle = igraph::layout_in_circle(graph_obj),
    igraph::layout_with_fr(graph_obj)
  )

  coords_dt <- as.data.table(coords)
  setnames(coords_dt, c("x", "y"))

  plot_nodes <- cbind(data.table::copy(node_dt), coords_dt)
  plot_nodes[, plot_label := ifelse(node_id %in% label_ids, gene_name, NA_character_)]
  plot_score_col <- if ("post_integrated" %in% names(plot_nodes)) {
    "post_integrated"
  } else if ("post_vulnerability" %in% names(plot_nodes)) {
    "post_vulnerability"
  } else {
    "post_norm"
  }
  plot_nodes[, plot_score := safe_numeric(get(plot_score_col))]

  edge_dt <- as.data.table(bundle$edges)
  edge_plot <- merge(
    edge_dt,
    plot_nodes[, .(from = node_id, x_from = x, y_from = y)],
    by = "from",
    all.x = TRUE
  )
  edge_plot <- merge(
    edge_plot,
    plot_nodes[, .(to = node_id, x_to = x, y_to = y)],
    by = "to",
    all.x = TRUE
  )

  edge_strength_col <- NULL
  if ("confidence_raw" %in% names(edge_plot)) {
    edge_strength_col <- "confidence_raw"
  } else if ("edge_reward_raw" %in% names(edge_plot)) {
    edge_strength_col <- "edge_reward_raw"
  } else if ("n_protein_edges" %in% names(edge_plot)) {
    edge_strength_col <- "n_protein_edges"
  }

  if (is.null(edge_strength_col)) {
    edge_plot[, edge_strength := 1]
  } else {
    edge_plot[, edge_strength := as.numeric(get(edge_strength_col))]
  }

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_plot,
      ggplot2::aes(
        x = x_from,
        y = y_from,
        xend = x_to,
        yend = y_to,
        alpha = edge_strength
      ),
      colour = "#3f4d63",
      linewidth = 0.45,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_nodes,
      ggplot2::aes(x = x, y = y, fill = plot_score),
      size = 6.2,
      shape = 21,
      colour = "#1f2937",
      stroke = 0.45
    ) +
    ggrepel::geom_label_repel(
      data = plot_nodes[!is.na(plot_label)],
      ggplot2::aes(x = x, y = y, label = plot_label),
      size = label_text_size,
      family = "Helvetica",
      fontface = "bold",
      label.size = 0.15,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.1, "lines"),
      fill = scales::alpha("white", 0.9),
      colour = "#111827",
      box.padding = label_box_padding,
      point.padding = label_point_padding,
      min.segment.length = 0,
      seed = 42,
      max.overlaps = Inf,
      force = label_force,
      max.time = 5,
      max.iter = 20000,
      segment.alpha = 0.35,
      segment.size = 0.2
    ) +
    ggplot2::scale_fill_gradient(
      low = "#d9e6f2",
      high = "#8a1538",
      name = "Post Integrated"
    ) +
    ggplot2::scale_alpha_continuous(
      range = c(0.12, 0.8),
      guide = "none"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void(base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, colour = "#111827", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "#4b5563"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.margin = grid::unit(c(18, 40, 18, 40), "pt")
    )
}

save_selected_subgraph_plot <- function(
  bundle,
  file_path,
  width = 12,
  height = 10,
  dpi = 300,
  layout = "fr",
  top_n_labels = Inf,
  title = "conseguiR Selected Gene Subgraph",
  subtitle = NULL
) {
  if (missing(file_path) || is.null(file_path) || !nzchar(file_path)) {
    stop("`file_path` must be provided to save the selected subgraph plot.")
  }

  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  plot_obj <- create_selected_subgraph_plot(
    bundle = bundle,
    layout = layout,
    top_n_labels = top_n_labels,
    title = title,
    subtitle = subtitle
  )

  ggplot2::ggsave(
    filename = file_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )

  invisible(file_path)
}

prepare_selected_subgraph_visualisation_bundle <- function(
  nodes = NULL,
  edges = NULL,
  summary = NULL,
  nodes_path = default_selected_subgraph_plot_config$nodes_path,
  edges_path = default_selected_subgraph_plot_config$edges_path,
  summary_path = default_selected_subgraph_plot_config$summary_path,
  sanitize_nonfinite = TRUE
) {
  if (is.null(nodes)) {
    nodes <- read_selected_subgraph_nodes(nodes_path)
  }
  if (is.null(edges)) {
    edges <- read_selected_subgraph_edges(edges_path)
  }
  if (is.null(summary)) {
    summary <- read_selected_subgraph_summary(summary_path)
  }

  summary_dt <- validate_selected_subgraph_summary(summary)
  plot_nodes <- prepare_subgraph_plot_nodes(nodes, sanitize_nonfinite = sanitize_nonfinite)
  plot_edges <- prepare_subgraph_plot_edges(edges, plot_nodes)
  graph <- build_selected_subgraph_igraph(plot_nodes, plot_edges)
  tidygraph_obj <- build_selected_subgraph_tidygraph(plot_nodes, plot_edges)
  cy_data <- prepare_cytoscape_data_frames(plot_nodes, plot_edges)

  list(
    nodes = plot_nodes,
    edges = plot_edges,
    summary = summary_dt,
    igraph = graph,
    tidygraph = tidygraph_obj,
    cytoscape_nodes = cy_data$nodes,
    cytoscape_edges = cy_data$edges
  )
}

save_selected_subgraph_visualisation_bundle <- function(
  bundle,
  output_prefix = default_selected_subgraph_plot_config$output_prefix,
  save_rds = TRUE,
  save_tables = TRUE
) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(save_rds)) {
    saveRDS(bundle, paste0(output_prefix, ".rds"))
  }

  if (isTRUE(save_tables)) {
    fwrite(bundle$nodes, paste0(output_prefix, "_nodes.tsv.gz"), sep = "\t")
    fwrite(bundle$edges, paste0(output_prefix, "_edges.tsv.gz"), sep = "\t")
    fwrite(bundle$summary, paste0(output_prefix, "_summary.tsv"), sep = "\t")
    fwrite(bundle$cytoscape_nodes, paste0(output_prefix, "_cytoscape_nodes.tsv.gz"), sep = "\t")
    fwrite(bundle$cytoscape_edges, paste0(output_prefix, "_cytoscape_edges.tsv.gz"), sep = "\t")
  }

  invisible(output_prefix)
}
