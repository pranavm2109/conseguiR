#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

default_config <- list(
  reg_elements_path = "data/raw/ENCODE/GRCh38-cCREs.bed",
  gene_links_zip_path = "data/raw/ENCODE/Human-Gene-Links.zip",
  gene_link_members = c(
    "V4-hg38.Gene-Links.3D-Chromatin.txt",
    "V4-hg38.Gene-Links.CRISPR.txt",
    "V4-hg38.Gene-Links.eQTLs.txt"
  ),
  gene_loc_path = "data/raw/NCBI38/NCBI38.gene.loc",
  output_prefix = "data/processed/gene_reg_graph_no_scores",
  min_link_value = 0,
  keep_self_loops = FALSE,
  directed = FALSE
)

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_encode_ccres <- function(path) {
  message("Reading ENCODE cCREs from: ", path)
  fread(path, header = FALSE, showProgress = FALSE)
}

read_encode_gene_links <- function(zip_path, member) {
  message("Reading ENCODE gene links from: ", basename(zip_path), "::", member)
  fread(
    cmd = sprintf("unzip -p %s %s", shQuote(zip_path), shQuote(member)),
    header = FALSE,
    showProgress = FALSE
  )
}

read_gene_loc_table <- function(path) {
  message("Reading gene locations from: ", path)
  fread(path, header = FALSE, showProgress = FALSE)
}

standardize_ccres <- function(reg_elements) {
  dt <- as.data.table(reg_elements)
  if (ncol(dt) < 6L) {
    stop(
      "ENCODE cCRE BED must have at least 6 columns: chrom, start, end, ",
      "dhs_id, ccre_id, ccre_class."
    )
  }

  unique(dt[, .(
    reg_chr = as.character(V1),
    reg_start = as.integer(V2) + 1L,
    reg_end = as.integer(V3),
    reg_accession = as.character(V4),
    reg_id = as.character(V5),
    reg_element_type = as.character(V6)
  )])[!is.na(reg_id) & reg_id != ""]
}

standardize_gene_loc_table <- function(gene_loc) {
  dt <- as.data.table(gene_loc)
  if (ncol(dt) < 6L) {
    stop("Gene location table must have at least 6 columns.")
  }

  dt <- dt[, .(
    gene_symbol = trimws(as.character(V6)),
    gene_chr = as.character(V2),
    gene_start = as.integer(V3),
    gene_end = as.integer(V4)
  )]
  dt <- dt[!is.na(gene_symbol) & gene_symbol != ""]
  dt[, .(
    gene_chr = first(gene_chr),
    gene_start = min(gene_start, na.rm = TRUE),
    gene_end = max(gene_end, na.rm = TRUE)
  ), by = gene_symbol]
}

source_from_member <- function(member) {
  if (grepl("3D-Chromatin", member, fixed = TRUE)) {
    return("3d_chromatin")
  }
  if (grepl("CRISPR", member, fixed = TRUE)) {
    return("crispr")
  }
  if (grepl("eQTL", member, fixed = TRUE)) {
    return("eqtl")
  }
  tolower(gsub("[^A-Za-z0-9]+", "_", member))
}

standardize_encode_links <- function(dt, source_name) {
  dt <- as.data.table(dt)

  if (identical(source_name, "3d_chromatin")) {
    if (ncol(dt) < 9L) {
      stop("ENCODE 3D chromatin links must have 9 columns.")
    }
    out <- dt[, .(
      reg_id = as.character(V1),
      ensembl_gene_id = as.character(V2),
      common_gene_name = trimws(as.character(V3)),
      gene_type = as.character(V4),
      evidence_type = source_name,
      assay_type = as.character(V5),
      experiment_id = as.character(V6),
      context_label = as.character(V7),
      effect_value = as.numeric(V8),
      p_value = suppressWarnings(as.numeric(V9))
    )]
  } else if (identical(source_name, "crispr")) {
    if (ncol(dt) < 10L) {
      stop("ENCODE CRISPR links must have 10 columns.")
    }
    out <- dt[, .(
      reg_id = as.character(V1),
      ensembl_gene_id = as.character(V2),
      common_gene_name = trimws(as.character(V3)),
      gene_type = as.character(V4),
      evidence_type = source_name,
      assay_type = as.character(V6),
      experiment_id = as.character(V7),
      context_label = as.character(V8),
      effect_value = as.numeric(V9),
      p_value = suppressWarnings(as.numeric(V10))
    )]
  } else if (identical(source_name, "eqtl")) {
    if (ncol(dt) < 9L) {
      stop("ENCODE eQTL links must have 9 columns.")
    }
    out <- dt[, .(
      reg_id = as.character(V1),
      ensembl_gene_id = as.character(V2),
      common_gene_name = trimws(as.character(V3)),
      gene_type = as.character(V4),
      evidence_type = source_name,
      assay_type = as.character(V6),
      experiment_id = as.character(V6),
      context_label = as.character(V7),
      effect_value = as.numeric(V8),
      p_value = suppressWarnings(as.numeric(V9))
    )]
  } else {
    stop("Unsupported ENCODE gene-link source: ", source_name)
  }

  out[, gene_id := fifelse(
    !is.na(common_gene_name) & common_gene_name != "",
    common_gene_name,
    ensembl_gene_id
  )]
  out[, score_component := fifelse(
    is.finite(p_value) & p_value > 0,
    -log10(pmax(p_value, 1e-300)),
    fifelse(is.finite(effect_value), abs(effect_value), 1)
  )]

  out[
    !is.na(reg_id) & reg_id != "" &
      !is.na(gene_id) & gene_id != "",
    .(
      reg_id,
      gene_id,
      ensembl_gene_id,
      gene_type,
      evidence_type,
      assay_type,
      experiment_id,
      context_label,
      effect_value,
      p_value,
      score_component
    )
  ]
}

read_all_encode_links <- function(zip_path, members) {
  rbindlist(
    lapply(members, function(member) {
      standardize_encode_links(
        read_encode_gene_links(zip_path, member),
        source_name = source_from_member(member)
      )
    }),
    fill = TRUE,
    use.names = TRUE
  )
}

collapse_encode_links <- function(links, reg_elements, gene_loc) {
  reg_dt <- as.data.table(reg_elements)
  gene_dt <- as.data.table(gene_loc)
  link_dt <- as.data.table(links)

  merged <- merge(link_dt, reg_dt, by = "reg_id", all.x = TRUE, sort = FALSE)
  merged <- merge(
    merged,
    gene_dt,
    by.x = "gene_id",
    by.y = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  )

  merged[
    ,
    .(
      ensembl_gene_id = first(ensembl_gene_id),
      gene_type = first(gene_type),
      link_value = uniqueN(evidence_type),
      link_score = sum(unique(score_component), na.rm = TRUE),
      link_method = paste(sort(unique(evidence_type)), collapse = "|"),
      evidence_count = .N,
      support_3d = any(evidence_type == "3d_chromatin"),
      support_crispr = any(evidence_type == "crispr"),
      support_eqtl = any(evidence_type == "eqtl"),
      assay_types = paste(sort(unique(na.omit(assay_type))), collapse = "|"),
      context_labels = paste(sort(unique(na.omit(context_label))), collapse = "|"),
      reg_chr = first(reg_chr),
      reg_start = first(reg_start),
      reg_end = first(reg_end),
      reg_accession = first(reg_accession),
      reg_element_type = first(reg_element_type),
      gene_chr = first(gene_chr),
      gene_start = first(gene_start),
      gene_end = first(gene_end)
    ),
    by = .(reg_id, gene_id)
  ]
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
    end = gene_end,
    ensembl_gene_id = ensembl_gene_id,
    gene_type = gene_type
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

  unique(
    rbindlist(list(gene_nodes, reg_nodes), fill = TRUE, use.names = TRUE),
    by = "node_id"
  )
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
    evidence_count,
    support_3d,
    support_crispr,
    support_eqtl,
    assay_types,
    context_labels,
    reg_chr,
    reg_start,
    reg_end,
    gene_chr,
    gene_start,
    gene_end
  )])
}

build_reg_target_labels <- function(edges) {
  dt <- as.data.table(edges)
  dt[
    ,
    .(label = paste(sort(unique(gene_id)), collapse = "|")),
    by = .(reg_elem_id = reg_id)
  ]
}

build_gene_reg_graph <- function(nodes, edges, directed = FALSE) {
  graph_from_data_frame(
    d = as.data.frame(edges),
    vertices = as.data.frame(nodes),
    directed = directed
  )
}

save_graph_outputs <- function(graph, nodes, edges, reg_target_labels, output_prefix) {
  ensure_parent_dir(output_prefix)

  saveRDS(graph, paste0(output_prefix, ".rds"))
  fwrite(nodes, paste0(output_prefix, "_nodes.tsv.gz"), sep = "\t")
  fwrite(edges, paste0(output_prefix, "_edges.tsv.gz"), sep = "\t")
  fwrite(
    reg_target_labels,
    file.path(dirname(output_prefix), "reg_target_labels.tsv.gz"),
    sep = "\t"
  )

  invisible(output_prefix)
}

prepare_gene_reg_graph <- function(config = default_config) {
  reg_elements_raw <- read_encode_ccres(config$reg_elements_path)
  gene_links_raw <- read_all_encode_links(
    zip_path = config$gene_links_zip_path,
    members = config$gene_link_members
  )
  gene_loc_raw <- read_gene_loc_table(config$gene_loc_path)

  reg_elements_std <- standardize_ccres(reg_elements_raw)
  gene_links_std <- collapse_encode_links(
    links = gene_links_raw,
    reg_elements = reg_elements_std,
    gene_loc = standardize_gene_loc_table(gene_loc_raw)
  )

  edges <- filter_gene_reg_links(
    edges = gene_links_std,
    min_link_value = config$min_link_value,
    keep_self_loops = config$keep_self_loops
  )

  nodes <- build_gene_reg_nodes(edges, reg_elements = reg_elements_std)
  edge_table <- build_gene_reg_edges(edges)
  reg_target_labels <- build_reg_target_labels(edges)
  graph <- build_gene_reg_graph(nodes, edge_table, directed = config$directed)

  save_graph_outputs(
    graph = graph,
    nodes = nodes,
    edges = edge_table,
    reg_target_labels = reg_target_labels,
    output_prefix = config$output_prefix
  )

  list(
    graph = graph,
    nodes = nodes,
    edges = edge_table,
    reg_target_labels = reg_target_labels,
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
