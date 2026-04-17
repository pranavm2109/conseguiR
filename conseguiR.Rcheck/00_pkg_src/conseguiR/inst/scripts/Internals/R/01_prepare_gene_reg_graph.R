#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})


default_config <- list(
  interactions_path = "data/raw/GeneHancer/gh_interactions_hg38_primary_assembly",
  reg_elements_path = "data/raw/GeneHancer/gh_reg_elements_hg38_primary_assembly",
  output_prefix = "data/processed/gene_reg_graph_no_scores",
  min_link_value = 0,
  keep_self_loops = FALSE,
  directed = FALSE
)

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_genehancer_interactions <- function(path) {
  message("Reading GeneHancer interactions from: ", path)
  fread(path)
}

read_genehancer_reg_elements <- function(path) {
  message("Reading GeneHancer regulatory elements from: ", path)
  fread(path)
}

standardize_genehancer_interactions <- function(interactions) {
  dt <- as.data.table(interactions)

  required_cols <- c(
    "geneHancerIdentifier",
    "geneName",
    "value",
    "score",
    "geneAssociationMethods",
    "geneHancerChrom",
    "geneHancerStart",
    "geneHancerEnd",
    "geneChrom",
    "geneStart",
    "geneEnd"
  )
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("GeneHancer interaction file is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  dt[, `:=`(
    reg_id = geneHancerIdentifier,
    gene_id = geneName,
    link_value = as.numeric(value),
    link_score = as.numeric(score),
    link_method = geneAssociationMethods,
    reg_chr = geneHancerChrom,
    reg_start = as.integer(geneHancerStart),
    reg_end = as.integer(geneHancerEnd),
    gene_chr = geneChrom,
    gene_start = as.integer(geneStart),
    gene_end = as.integer(geneEnd)
  )]

  dt <- unique(
    dt[, .(
      reg_id,
      gene_id,
      link_value,
      link_score,
      link_method,
      reg_chr,
      reg_start,
      reg_end,
      gene_chr,
      gene_start,
      gene_end
    )]
  )

  dt[!is.na(reg_id) & reg_id != "" & !is.na(gene_id) & gene_id != ""]
}

standardize_reg_elements <- function(reg_elements) {
  dt <- as.data.table(reg_elements)

  required_cols <- c("name", "#chrom", "chromStart", "chromEnd", "score", "elementType", "eliteness")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("GeneHancer reg-elements file is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  unique(dt[, .(
    reg_id = name,
    reg_chr = get("#chrom"),
    reg_start = as.integer(chromStart),
    reg_end = as.integer(chromEnd),
    reg_score = as.numeric(score),
    reg_element_type = elementType,
    reg_eliteness = eliteness
  )])
}

filter_gene_reg_links <- function(edges, min_link_value = 0, keep_self_loops = FALSE) {
  dt <- copy(as.data.table(edges))
  dt <- dt[is.na(link_value) | link_value >= min_link_value]

  if (!keep_self_loops) {
    dt <- dt[reg_id != gene_id]
  }

  unique(dt)
}

build_gene_reg_nodes <- function(edges, reg_elements = NULL) {
  edge_dt <- as.data.table(edges)

  gene_nodes <- unique(edge_dt[, .(
    name = gene_id,
    node_id = gene_id,
    node_type = "gene",
    chr = gene_chr,
    start = gene_start,
    end = gene_end
  )])

  reg_nodes_from_edges <- unique(edge_dt[, .(
    name = reg_id,
    node_id = reg_id,
    node_type = "reg",
    chr = reg_chr,
    start = reg_start,
    end = reg_end
  )])

  if (!is.null(reg_elements)) {
    reg_dt <- as.data.table(reg_elements)
    reg_nodes <- merge(
      reg_nodes_from_edges,
      reg_dt,
      by.x = "node_id",
      by.y = "reg_id",
      all.x = TRUE,
      sort = FALSE
    )
    setnames(reg_nodes, "name", "node_name")
    reg_nodes[, name := node_id]
    reg_nodes[, node_name := NULL]
  } else {
    reg_nodes <- reg_nodes_from_edges
  }

  rbindlist(list(gene_nodes, reg_nodes), fill = TRUE, use.names = TRUE)
}

build_gene_reg_edges <- function(edges) {
  dt <- copy(as.data.table(edges))

  dt[, `:=`(
    from = reg_id,
    to = gene_id,
    weight = fifelse(is.na(link_value), NA_real_, 1 / (1 + link_value)),
    confidence = link_value
  )]

  unique(dt[, .(
    from,
    to,
    weight,
    confidence,
    link_score,
    link_method,
    reg_chr,
    reg_start,
    reg_end,
    gene_chr,
    gene_start,
    gene_end
  )])
}

build_gene_reg_graph <- function(nodes, edges, directed = FALSE) {
  graph_from_data_frame(
    d = as.data.frame(edges),
    vertices = as.data.frame(nodes),
    directed = directed
  )
}

save_graph_outputs <- function(graph, nodes, edges, output_prefix) {
  ensure_parent_dir(output_prefix)

  saveRDS(graph, paste0(output_prefix, ".rds"))
  fwrite(nodes, paste0(output_prefix, "_nodes.tsv.gz"), sep = "\t")
  fwrite(edges, paste0(output_prefix, "_edges.tsv.gz"), sep = "\t")

  invisible(output_prefix)
}

prepare_gene_reg_graph <- function(config = default_config) {
  interactions_raw <- read_genehancer_interactions(config$interactions_path)
  reg_elements_raw <- read_genehancer_reg_elements(config$reg_elements_path)

  interactions_std <- standardize_genehancer_interactions(interactions_raw)
  reg_elements_std <- standardize_reg_elements(reg_elements_raw)

  edges <- filter_gene_reg_links(
    edges = interactions_std,
    min_link_value = config$min_link_value,
    keep_self_loops = config$keep_self_loops
  )

  nodes <- build_gene_reg_nodes(edges, reg_elements = reg_elements_std)
  edge_table <- build_gene_reg_edges(edges)
  graph <- build_gene_reg_graph(nodes, edge_table, directed = config$directed)

  save_graph_outputs(
    graph = graph,
    nodes = nodes,
    edges = edge_table,
    output_prefix = config$output_prefix
  )

  list(
    graph = graph,
    nodes = nodes,
    edges = edge_table,
    config = config
  )
}

main <- function() {
  result <- prepare_gene_reg_graph()

  message("Gene-reg graph preparation complete.")
  message("Nodes: ", nrow(result$nodes))
  message("Edges: ", nrow(result$edges))
  invisible(result)
}

if (sys.nframe() == 0) {
  main()
}
