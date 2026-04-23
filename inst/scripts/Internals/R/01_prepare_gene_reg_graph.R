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
  evidence_alpha = 0.6,
  weight_3d = 1.0,
  weight_crispr = 1.25,
  weight_eqtl = 1.0,
  support_bonus = 0.15,
  min_link_value = 0,
  keep_self_loops = FALSE,
  directed = FALSE
)

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_compact_table_xz <- function(dt, path) {
  ensure_parent_dir(path)
  con <- xzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    as.data.frame(dt),
    file = con,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  invisible(path)
}

read_encode_ccres <- function(path) {
  message("Reading ENCODE cCREs from: ", path)
  data.table::fread(path, header = FALSE, showProgress = FALSE)
}

read_encode_gene_links <- function(zip_path, member) {
  message("Reading ENCODE gene links from: ", basename(zip_path), "::", member)
  extract_dir <- tempfile("conseguiR_encode_links_")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip_path, files = member, exdir = extract_dir)
  extracted_path <- file.path(extract_dir, member)
  if (!file.exists(extracted_path)) {
    stop("Failed to extract ENCODE gene-link member: ", member)
  }
  data.table::as.data.table(
    utils::read.delim(
      extracted_path,
      header = FALSE,
      sep = "\t",
      fill = TRUE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}

read_gene_loc_table <- function(path) {
  message("Reading gene locations from: ", path)
  data.table::fread(path, header = FALSE, showProgress = FALSE)
}

standardize_ccres <- function(reg_elements) {
  dt <- as.data.table(reg_elements)
  if (ncol(dt) < 6L) {
    stop(
      "ENCODE cCRE BED must have at least 6 columns: chrom, start, end, ",
      "dhs_id, ccre_id, ccre_class."
    )
  }

  unique(dt[, list(
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

  dt <- dt[, list(
    gene_symbol = trimws(as.character(V6)),
    gene_chr = as.character(V2),
    gene_start = as.integer(V3),
    gene_end = as.integer(V4)
  )]
  dt <- dt[!is.na(gene_symbol) & gene_symbol != ""]
  dt[, list(
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

percentile_rank <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) {
    return(out)
  }
  if (sum(ok) == 1L) {
    out[ok] <- 1
    return(out)
  }
  ranks <- data.table::frank(x[ok], ties.method = "average", na.last = "keep")
  out[ok] <- (ranks - 1) / (sum(ok) - 1)
  out
}

max_finite_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  max(x)
}

first_non_missing_character <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) {
    return(NA_character_)
  }
  as.character(x[[1]])
}

standardize_encode_links <- function(dt, source_name, evidence_alpha = 0.6) {
  dt <- as.data.table(dt)

  if (identical(source_name, "3d_chromatin")) {
    if (ncol(dt) < 9L) {
      stop("ENCODE 3D chromatin links must have 9 columns.")
    }
    out <- dt[, list(
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
    out <- dt[, list(
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
    out <- dt[, list(
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
  out[, metric_value := as.numeric(effect_value)]
  out[, metric_magnitude := abs(metric_value)]
  out[, significance_value := fifelse(
    is.finite(p_value) & p_value > 0,
    -log10(pmax(p_value, 1e-300)),
    NA_real_
  )]
  out[, source_magnitude_rank := percentile_rank(metric_magnitude)]
  out[, source_significance_rank := percentile_rank(significance_value)]
  out[, source_evidence := data.table::fifelse(
    is.finite(source_significance_rank) & is.finite(source_magnitude_rank),
    evidence_alpha * source_significance_rank +
      (1 - evidence_alpha) * source_magnitude_rank,
    data.table::fifelse(
      is.finite(source_significance_rank),
      source_significance_rank,
      source_magnitude_rank
    )
  )]

  out[
    !is.na(reg_id) & reg_id != "" &
      !is.na(gene_id) & gene_id != "",
    list(
      reg_id,
      gene_id,
      ensembl_gene_id,
      gene_type,
      evidence_type,
      assay_type,
      experiment_id,
      context_label,
      metric_value,
      metric_magnitude,
      p_value,
      significance_value,
      source_significance_rank,
      source_magnitude_rank,
      source_evidence
    )
  ]
}

collapse_encode_links <- function(links, reg_elements, gene_loc, source_name) {
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

  collapsed <- merged[
    ,
    list(
      ensembl_gene_id = first(ensembl_gene_id),
      gene_type = first(gene_type),
      source_row_count = .N,
      metric_value = metric_value[which.max(data.table::fifelse(
        is.finite(source_evidence),
        source_evidence,
        -Inf
      ))],
      p_value = p_value[which.max(data.table::fifelse(
        is.finite(source_evidence),
        source_evidence,
        -Inf
      ))],
      source_significance = max_finite_or_na(source_significance_rank),
      source_magnitude = max_finite_or_na(source_magnitude_rank),
      source_evidence = max_finite_or_na(source_evidence),
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
    by = list(reg_id, gene_id)
  ]

  for (col in c("source_significance", "source_magnitude", "source_evidence")) {
    collapsed[!is.finite(get(col)), (col) := NA_real_]
  }

  metric_col <- switch(
    source_name,
    "3d_chromatin" = "score_3d",
    "crispr" = "effect_crispr",
    "eqtl" = "slope_eqtl",
    "metric_value"
  )
  p_col <- switch(
    source_name,
    "3d_chromatin" = "p_3d",
    "crispr" = "p_crispr",
    "eqtl" = "p_eqtl",
    "p_value"
  )
  sig_col <- paste0("significance_", source_name)
  mag_col <- paste0("magnitude_", source_name)
  evidence_col <- paste0("evidence_", source_name)
  assay_col <- paste0("assay_types_", source_name)
  context_col <- paste0("context_labels_", source_name)
  rows_col <- paste0("rows_", source_name)
  support_col <- switch(
    source_name,
    "3d_chromatin" = "support_3d",
    "crispr" = "support_crispr",
    "eqtl" = "support_eqtl",
    paste0("support_", source_name)
  )

  data.table::setnames(collapsed, "metric_value", metric_col)
  data.table::setnames(collapsed, "p_value", p_col)
  data.table::setnames(collapsed, "source_significance", sig_col)
  data.table::setnames(collapsed, "source_magnitude", mag_col)
  data.table::setnames(collapsed, "source_evidence", evidence_col)
  data.table::setnames(collapsed, "assay_types", assay_col)
  data.table::setnames(collapsed, "context_labels", context_col)
  data.table::setnames(collapsed, "source_row_count", rows_col)
  collapsed[, (support_col) := TRUE]
  collapsed
}

process_encode_link_member <- function(
  member,
  zip_path,
  reg_elements,
  gene_loc,
  evidence_alpha = 0.6
) {
  source_name <- source_from_member(member)
  raw_links <- read_encode_gene_links(zip_path, member)
  if (is.null(raw_links) || nrow(raw_links) == 0L) {
    return(NULL)
  }
  collapse_encode_links(
    links = standardize_encode_links(
      raw_links,
      source_name = source_name,
      evidence_alpha = evidence_alpha
    ),
    reg_elements = reg_elements,
    gene_loc = gene_loc,
    source_name = source_name
  )
}

collapse_all_encode_links <- function(
  zip_path,
  members,
  reg_elements,
  gene_loc,
  evidence_alpha = 0.6,
  weight_3d = 1,
  weight_crispr = 1.25,
  weight_eqtl = 1,
  support_bonus = 0.15
) {
  partials <- lapply(
    members,
    process_encode_link_member,
    zip_path = zip_path,
    reg_elements = reg_elements,
    gene_loc = gene_loc,
    evidence_alpha = evidence_alpha
  )
  partials <- Filter(Negate(is.null), partials)
  if (length(partials) == 0L) {
    return(data.table::data.table())
  }
  combined <- Reduce(function(x, y) {
    merge(x, y, by = c("reg_id", "gene_id"), all = TRUE, sort = FALSE)
  }, partials)
  combined <- data.table::as.data.table(combined)

  required_columns <- c(
    "ensembl_gene_id", "gene_type",
    "score_3d", "p_3d", "significance_3d_chromatin",
    "magnitude_3d_chromatin", "evidence_3d_chromatin",
    "effect_crispr", "p_crispr", "significance_crispr",
    "magnitude_crispr", "evidence_crispr",
    "slope_eqtl", "p_eqtl", "significance_eqtl",
    "magnitude_eqtl", "evidence_eqtl",
    "rows_3d_chromatin", "rows_crispr", "rows_eqtl",
    "support_3d", "support_crispr", "support_eqtl",
    "assay_types_3d_chromatin", "assay_types_crispr", "assay_types_eqtl",
    "context_labels_3d_chromatin", "context_labels_crispr",
    "context_labels_eqtl", "reg_chr", "reg_start", "reg_end",
    "reg_accession", "reg_element_type", "gene_chr", "gene_start", "gene_end"
  )
  missing_columns <- setdiff(required_columns, names(combined))
  for (col in missing_columns) {
    combined[, (col) := NA]
  }

  combined[
    ,
    list(
      ensembl_gene_id = first_non_missing_character(ensembl_gene_id),
      gene_type = first_non_missing_character(gene_type),
      score_3d = first(score_3d),
      p_3d = first(p_3d),
      significance_3d_chromatin = first(significance_3d_chromatin),
      magnitude_3d_chromatin = first(magnitude_3d_chromatin),
      evidence_3d_chromatin = first(evidence_3d_chromatin),
      effect_crispr = first(effect_crispr),
      p_crispr = first(p_crispr),
      significance_crispr = first(significance_crispr),
      magnitude_crispr = first(magnitude_crispr),
      evidence_crispr = first(evidence_crispr),
      slope_eqtl = first(slope_eqtl),
      p_eqtl = first(p_eqtl),
      significance_eqtl = first(significance_eqtl),
      magnitude_eqtl = first(magnitude_eqtl),
      evidence_eqtl = first(evidence_eqtl),
      rows_3d_chromatin = first(rows_3d_chromatin),
      rows_crispr = first(rows_crispr),
      rows_eqtl = first(rows_eqtl),
      support_3d = isTRUE(first(support_3d)),
      support_crispr = isTRUE(first(support_crispr)),
      support_eqtl = isTRUE(first(support_eqtl)),
      assay_types_3d_chromatin = first(assay_types_3d_chromatin),
      assay_types_crispr = first(assay_types_crispr),
      assay_types_eqtl = first(assay_types_eqtl),
      context_labels_3d_chromatin = first(context_labels_3d_chromatin),
      context_labels_crispr = first(context_labels_crispr),
      context_labels_eqtl = first(context_labels_eqtl),
      reg_chr = first(reg_chr),
      reg_start = first(reg_start),
      reg_end = first(reg_end),
      reg_accession = first(reg_accession),
      reg_element_type = first(reg_element_type),
      gene_chr = first(gene_chr),
      gene_start = first(gene_start),
      gene_end = first(gene_end)
    ),
    by = list(reg_id, gene_id)
  ][
    ,
    `:=`(
      support_count = as.integer(support_3d) +
        as.integer(support_crispr) +
        as.integer(support_eqtl),
      combined_edge_score =
        weight_3d * data.table::fcoalesce(evidence_3d_chromatin, 0) +
        weight_crispr * data.table::fcoalesce(evidence_crispr, 0) +
        weight_eqtl * data.table::fcoalesce(evidence_eqtl, 0) +
        support_bonus * pmax(
          as.integer(support_3d) +
            as.integer(support_crispr) +
            as.integer(support_eqtl) - 1L,
          0L
        )
    )
  ][
    ,
    `:=`(
      link_value = combined_edge_score,
      link_score = combined_edge_score,
      evidence_count = data.table::fcoalesce(rows_3d_chromatin, 0L) +
        data.table::fcoalesce(rows_crispr, 0L) +
        data.table::fcoalesce(rows_eqtl, 0L)
    )
  ][
    ,
    link_method := data.table::fifelse(
      support_3d & support_crispr & support_eqtl,
      "3d_chromatin|crispr|eqtl",
      data.table::fifelse(
        support_3d & support_crispr,
        "3d_chromatin|crispr",
        data.table::fifelse(
          support_3d & support_eqtl,
          "3d_chromatin|eqtl",
          data.table::fifelse(
            support_crispr & support_eqtl,
            "crispr|eqtl",
            data.table::fifelse(
              support_3d,
              "3d_chromatin",
              data.table::fifelse(
                support_crispr,
                "crispr",
                data.table::fifelse(
                  support_eqtl,
                  "eqtl",
                  NA_character_
                )
              )
            )
          )
        )
      )
    )]
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

  gene_nodes <- unique(edge_dt[, list(
    name = gene_id,
    node_id = gene_id,
    node_type = "gene",
    chr = gene_chr,
    start = gene_start,
    end = gene_end,
    ensembl_gene_id = ensembl_gene_id,
    gene_type = gene_type
  )])

  reg_nodes_from_edges <- unique(edge_dt[, list(
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
  )[, node_index := seq_len(.N)]
}

build_gene_reg_edges <- function(edges) {
  dt <- copy(as.data.table(edges))

  dt[, `:=`(
    from = reg_id,
    to = gene_id,
    confidence = combined_edge_score,
    weight = fifelse(
      is.na(combined_edge_score),
      NA_real_,
      1 / (1 + combined_edge_score)
    )
  )]
  dt[, link_method := data.table::fifelse(
    support_3d & support_crispr & support_eqtl,
    "3d_chromatin|crispr|eqtl",
    data.table::fifelse(
      support_3d & support_crispr,
      "3d_chromatin|crispr",
      data.table::fifelse(
        support_3d & support_eqtl,
        "3d_chromatin|eqtl",
        data.table::fifelse(
          support_crispr & support_eqtl,
          "crispr|eqtl",
          data.table::fifelse(
            support_3d,
            "3d_chromatin",
            data.table::fifelse(
              support_crispr,
              "crispr",
              data.table::fifelse(
                support_eqtl,
                "eqtl",
                NA_character_
              )
            )
          )
        )
      )
    )
  )]

  unique(dt[, list(
    from,
    to,
    weight,
    confidence,
    combined_edge_score,
    link_score,
    link_method,
    support_count,
    evidence_count,
    score_3d,
    p_3d,
    significance_3d_chromatin,
    magnitude_3d_chromatin,
    evidence_3d_chromatin,
    effect_crispr,
    p_crispr,
    significance_crispr,
    magnitude_crispr,
    evidence_crispr,
    slope_eqtl,
    p_eqtl,
    significance_eqtl,
    magnitude_eqtl,
    evidence_eqtl,
    support_3d,
    support_crispr,
    support_eqtl,
    rows_3d_chromatin,
    rows_crispr,
    rows_eqtl,
    assay_types_3d_chromatin,
    assay_types_crispr,
    assay_types_eqtl,
    context_labels_3d_chromatin,
    context_labels_crispr,
    context_labels_eqtl,
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
    list(label = paste(sort(unique(gene_id)), collapse = "|")),
    by = list(reg_elem_id = reg_id)
  ]
}

compact_gene_reg_nodes <- function(nodes) {
  dt <- copy(as.data.table(nodes))
  if (!"node_index" %in% names(dt)) {
    dt[, node_index := seq_len(.N)]
  }
  dt
}

compact_gene_reg_edges <- function(edges, nodes) {
  edge_dt <- copy(as.data.table(edges))
  node_dt <- compact_gene_reg_nodes(nodes)[, .(node_index, node_id)]
  from_lookup <- node_dt[, .(from = node_id, from_idx = node_index)]
  to_lookup <- node_dt[, .(to = node_id, to_idx = node_index)]

  edge_dt <- merge(edge_dt, from_lookup, by = "from", all.x = TRUE, sort = FALSE)
  edge_dt <- merge(edge_dt, to_lookup, by = "to", all.x = TRUE, sort = FALSE)
  edge_dt[, .(
    from_idx = as.integer(from_idx),
    to_idx = as.integer(to_idx),
    confidence = as.numeric(confidence),
    link_method = as.character(link_method)
  )]
}

build_gene_reg_graph <- function(nodes, edges, directed = FALSE) {
  vertex_df <- as.data.frame(nodes[, !"node_index"])
  graph_from_data_frame(
    d = as.data.frame(edges),
    vertices = vertex_df,
    directed = directed
  )
}

save_graph_outputs <- function(graph, nodes, edges, reg_target_labels, output_prefix) {
  ensure_parent_dir(output_prefix)

  saveRDS(graph, paste0(output_prefix, ".rds"))
  fwrite(nodes[, !"node_index"], paste0(output_prefix, "_nodes.tsv.gz"), sep = "\t")
  fwrite(edges, paste0(output_prefix, "_edges.tsv.gz"), sep = "\t")
  write_compact_table_xz(compact_gene_reg_nodes(nodes), paste0(output_prefix, "_nodes_compact.tsv.xz"))
  write_compact_table_xz(compact_gene_reg_edges(edges, nodes), paste0(output_prefix, "_edges_compact.tsv.xz"))
  fwrite(
    reg_target_labels,
    file.path(dirname(output_prefix), "reg_target_labels.tsv.gz"),
    sep = "\t"
  )

  invisible(output_prefix)
}

prepare_gene_reg_graph <- function(config = default_config) {
  config <- utils::modifyList(default_config, config)

  reg_elements_raw <- read_encode_ccres(config$reg_elements_path)
  gene_loc_raw <- read_gene_loc_table(config$gene_loc_path)

  reg_elements_std <- standardize_ccres(reg_elements_raw)
  gene_loc_std <- standardize_gene_loc_table(gene_loc_raw)
  gene_links_std <- collapse_all_encode_links(
    zip_path = config$gene_links_zip_path,
    members = config$gene_link_members,
    reg_elements = reg_elements_std,
    gene_loc = gene_loc_std,
    evidence_alpha = config$evidence_alpha,
    weight_3d = config$weight_3d,
    weight_crispr = config$weight_crispr,
    weight_eqtl = config$weight_eqtl,
    support_bonus = config$support_bonus
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
