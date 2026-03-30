#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

default_config <- list(
  protein_links_path = "data/raw/STRING/9606.protein.links.v12.0.txt",
  protein_info_path = "data/raw/STRING/9606.protein.info.v12.0.txt",
  output_prefix = "data/processed/gene_gene_graph",
  min_combined_score = 400,
  directed = FALSE,
  collapse_to_gene_level = TRUE
)

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_string_links <- function(path) {
  message("Reading STRING links from: ", path)
  fread(path)
}

read_string_protein_info <- function(path) {
  message("Reading STRING protein metadata from: ", path)
  fread(path)
}

standardize_string_links <- function(links) {
  dt <- as.data.table(links)

  required_cols <- c("protein1", "protein2", "combined_score")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("STRING links file is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  unique(dt[, .(
    protein1 = protein1,
    protein2 = protein2,
    combined_score = as.numeric(combined_score)
  )])
}

standardize_string_protein_info <- function(protein_info) {
  dt <- as.data.table(protein_info)

  required_cols <- c("#string_protein_id", "preferred_name")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("STRING protein info file is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  unique(dt[, .(
    protein_id = get("#string_protein_id"),
    gene_symbol = preferred_name,
    protein_size = if ("protein_size" %in% names(dt)) as.integer(protein_size) else NA_integer_,
    annotation = if ("annotation" %in% names(dt)) annotation else NA_character_
  )])
}

filter_string_links <- function(links, min_combined_score = 400) {
  dt <- as.data.table(links)
  dt[!is.na(combined_score) & combined_score >= min_combined_score]
}

map_links_to_genes <- function(links, protein_info) {
  info_dt <- as.data.table(protein_info)
  link_dt <- as.data.table(links)

  mapped <- merge(link_dt, info_dt, by.x = "protein1", by.y = "protein_id", all.x = TRUE, sort = FALSE)
  setnames(mapped, c("gene_symbol", "protein_size", "annotation"), c("gene1", "protein1_size", "protein1_annotation"))

  mapped <- merge(mapped, info_dt, by.x = "protein2", by.y = "protein_id", all.x = TRUE, sort = FALSE)
  setnames(mapped, c("gene_symbol", "protein_size", "annotation"), c("gene2", "protein2_size", "protein2_annotation"))

  mapped
}

collapse_string_to_gene_level <- function(mapped_links) {
  dt <- copy(as.data.table(mapped_links))
  dt <- dt[!is.na(gene1) & gene1 != "" & !is.na(gene2) & gene2 != ""]

  dt[, c("gene_a", "gene_b") := .(
    pmin(gene1, gene2),
    pmax(gene1, gene2)
  )]

  dt <- dt[gene_a != gene_b]

  dt[, .(
    combined_score = max(combined_score, na.rm = TRUE),
    n_protein_edges = .N
  ), by = .(gene_a, gene_b)]
}

build_gene_gene_nodes <- function(mapped_links) {
  dt <- as.data.table(mapped_links)
  gene_ids <- sort(unique(c(dt$gene1, dt$gene2)))
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]

  data.table(
    name = gene_ids,
    node_id = gene_ids,
    node_type = "gene"
  )
}

build_gene_gene_edges <- function(gene_level_links) {
  dt <- copy(as.data.table(gene_level_links))

  dt[, `:=`(
    from = gene_a,
    to = gene_b,
    confidence = combined_score,
    weight = 1 / (1 + combined_score)
  )]

  unique(dt[, .(
    from,
    to,
    confidence,
    weight,
    n_protein_edges
  )])
}

build_gene_gene_graph <- function(nodes, edges, directed = FALSE) {
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

prepare_gene_gene_graph <- function(config = default_config) {
  links_raw <- read_string_links(config$protein_links_path)
  protein_info_raw <- read_string_protein_info(config$protein_info_path)

  links_std <- standardize_string_links(links_raw)
  protein_info_std <- standardize_string_protein_info(protein_info_raw)
  links_filtered <- filter_string_links(links_std, min_combined_score = config$min_combined_score)

  mapped_links <- map_links_to_genes(links_filtered, protein_info_std)

  if (isTRUE(config$collapse_to_gene_level)) {
    gene_level_links <- collapse_string_to_gene_level(mapped_links)
  } else {
    stop("Only collapse_to_gene_level = TRUE is currently scaffolded.")
  }

  nodes <- build_gene_gene_nodes(mapped_links)
  edges <- build_gene_gene_edges(gene_level_links)
  graph <- build_gene_gene_graph(nodes, edges, directed = config$directed)

  save_graph_outputs(
    graph = graph,
    nodes = nodes,
    edges = edges,
    output_prefix = config$output_prefix
  )

  list(
    graph = graph,
    nodes = nodes,
    edges = edges,
    config = config
  )
}

main <- function() {
  result <- prepare_gene_gene_graph()

  message("Gene-gene graph preparation complete.")
  message("Nodes: ", nrow(result$nodes))
  message("Edges: ", nrow(result$edges))
  invisible(result)
}

if (sys.nframe() == 0) {
  main()
}
