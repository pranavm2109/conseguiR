#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

conseguiR_plot_runtime_file <- function(relpath) {
  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(candidate)
  }

  pkg_candidates <- unique(c(
    relpath,
    sub("^inst/", "", relpath)
  ))
  for (pkg_relpath in pkg_candidates) {
    pkg_path <- system.file(pkg_relpath, package = "conseguiR")
    if (nzchar(pkg_path) && file.exists(pkg_path)) {
      return(pkg_path)
    }
  }

  stop("Could not locate required plotting runtime file: ", relpath)
}

default_stage_plot_config <- list(
  score_output_prefix = "data/processed/conseguiR_score_plot"
)

.conseguiR_plot_cache <- new.env(parent = emptyenv())

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

read_backend_gene_label_map <- function() {
  if (exists("gene_label_map", envir = .conseguiR_plot_cache, inherits = FALSE)) {
    return(data.table::copy(get("gene_label_map", envir = .conseguiR_plot_cache, inherits = FALSE)))
  }

  loc_path <- conseguiR_plot_runtime_file("inst/extdata/backend/NCBI38.gene.loc")
  dt <- data.table::fread(loc_path, header = FALSE, showProgress = FALSE)
  if (ncol(dt) < 6L) {
    stop("Backend gene location file does not contain the expected label column.")
  }

  out <- unique(dt[, .(
    feature_id = as.character(V1),
    label = as.character(V6)
  )], by = "feature_id")
  out <- out[!is.na(feature_id) & feature_id != "" & !is.na(label) & label != ""]
  assign("gene_label_map", out, envir = .conseguiR_plot_cache)
  data.table::copy(out)
}

read_backend_reg_label_map <- function() {
  if (exists("reg_label_map", envir = .conseguiR_plot_cache, inherits = FALSE)) {
    return(data.table::copy(get("reg_label_map", envir = .conseguiR_plot_cache, inherits = FALSE)))
  }

  mapping_path <- conseguiR_plot_runtime_file("inst/extdata/backend/genehancer_reg_target_labels.tsv.gz")
  if (!file.exists(mapping_path)) {
    return(data.table(reg_elem_id = character(), label = character()))
  }

  dt <- data.table::as.data.table(data.table::fread(mapping_path, showProgress = FALSE))
  out <- unique(dt[!is.na(reg_elem_id) & reg_elem_id != "" & !is.na(label) & label != ""], by = "reg_elem_id")
  assign("reg_label_map", out, envir = .conseguiR_plot_cache)
  data.table::copy(out)
}

resolve_score_feature_column <- function(dt) {
  candidates <- c("gene_id", "feature_id", "reg_elem_id", "gene_name")
  match <- intersect(candidates, names(dt))
  if (length(match) == 0L) {
    stop("Could not infer a feature column for score plotting.")
  }
  match[[1]]
}

infer_feature_context <- function(feature_column, which = NULL, dt = NULL) {
  if (identical(feature_column, "reg_elem_id")) {
    return("reg")
  }
  if (identical(feature_column, "gene_id") || identical(feature_column, "gene_name")) {
    return("gene")
  }
  if (!is.null(which) && grepl("reg", which, ignore.case = TRUE)) {
    return("reg")
  }
  if (!is.null(which) && grepl("gene", which, ignore.case = TRUE)) {
    return("gene")
  }
  if (!is.null(dt) && "feature_id" %in% names(dt)) {
    feature_values <- unique(stats::na.omit(as.character(dt$feature_id)))
    feature_values <- feature_values[nzchar(feature_values)]
    if (length(feature_values) > 0L && any(grepl("^GH", feature_values))) {
      return("reg")
    }
  }
  "gene"
}

resolve_score_label_map <- function(feature_column, which = NULL, dt = NULL) {
  feature_context <- infer_feature_context(
    feature_column = feature_column,
    which = which,
    dt = dt
  )

  switch(
    feature_context,
    reg = read_backend_reg_label_map(),
    gene = read_backend_gene_label_map(),
    NULL
  )
}

prepare_score_plot_table <- function(
  table,
  which = NULL,
  feature_column = NULL,
  z_column = "zstat",
  p_value_column = NULL,
  label_features = NULL
) {
  dt <- data.table::as.data.table(data.table::copy(table))
  feature_column <- feature_column %||% resolve_score_feature_column(dt)

  if (!feature_column %in% names(dt)) {
    stop("Feature column not found in score table: ", feature_column)
  }
  if (!z_column %in% names(dt)) {
    stop("Z-score column not found in score table: ", z_column)
  }

  dt[, feature_id_plot := as.character(get(feature_column))]
  dt[, z_raw := as.numeric(get(z_column))]
  dt <- dt[!is.na(feature_id_plot) & feature_id_plot != "" & !is.na(z_raw)]

  finite_z <- dt[is.finite(z_raw), z_raw]
  z_cap <- if (length(finite_z) > 0L) max(abs(finite_z), na.rm = TRUE) + 0.5 else 10
  dt[, z_plot := z_raw]
  dt[is.infinite(z_plot) & z_plot > 0, z_plot := z_cap]
  dt[is.infinite(z_plot) & z_plot < 0, z_plot := -z_cap]

  label_map <- resolve_score_label_map(
    feature_column = feature_column,
    which = which,
    dt = dt
  )
  if (!is.null(label_map) && nrow(label_map) > 0L) {
    merge_cols <- names(label_map)[1]
    dt <- merge(
      dt,
      label_map,
      by.x = "feature_id_plot",
      by.y = merge_cols,
      all.x = TRUE
    )
    dt[, feature_label := data.table::fifelse(!is.na(label) & label != "", label, feature_id_plot)]
    dt[, label := NULL]
  } else {
    dt[, feature_label := feature_id_plot]
  }

  if (!is.null(p_value_column) && p_value_column %in% names(dt)) {
    dt[, p_value_plot := as.numeric(get(p_value_column))]
  } else if ("p_value" %in% names(dt)) {
    dt[, p_value_plot := as.numeric(p_value)]
  } else {
    dt[, p_value_plot := NA_real_]
  }

  if (all(!is.finite(dt$p_value_plot))) {
    dt[, p_value_plot := 2 * stats::pnorm(-abs(z_plot))]
  }

  dt[p_value_plot <= 0 | !is.finite(p_value_plot), p_value_plot := 1e-300]

  dt[, neglog10_p := ifelse(
    is.finite(p_value_plot) & p_value_plot > 0,
    -log10(pmax(p_value_plot, 1e-300)),
    NA_real_
  )]

  data.table::setorder(dt, -z_plot)
  dt[, rank_plot := seq_len(.N)]
  requested_labels <- unique(as.character(label_features %||% character()))
  requested_labels <- requested_labels[nzchar(requested_labels)]

  dt[, feature_label_tokens := strsplit(feature_label, "|", fixed = TRUE)]
  dt[, matched_label := vapply(
    feature_label_tokens,
    function(labels) {
      matched <- intersect(labels, requested_labels)
      if (length(matched) == 0L) "" else matched[[1]]
    },
    character(1)
  )]
  dt[, highlighted := feature_id_plot %in% requested_labels |
       feature_label %in% requested_labels |
       matched_label != ""]
  dt[, should_label := highlighted]
  dt[matched_label != "",
     label_rank := data.table::frank(-z_plot, ties.method = "first"),
     by = matched_label]
  dt[matched_label != "" & label_rank > 1L, should_label := FALSE]
  dt[, feature_label_tokens := NULL]
  dt[, label_rank := NULL]
  dt
}

infer_plot_mode <- function(dt, test_tail = c("auto", "one_tailed", "two_tailed")) {
  test_tail <- match.arg(test_tail)
  if (identical(test_tail, "one_tailed")) {
    return("rank")
  }
  if (identical(test_tail, "two_tailed")) {
    return("volcano")
  }

  if ("p_value_plot" %in% names(dt) && any(is.finite(dt$p_value_plot))) {
    "volcano"
  } else {
    "rank"
  }
}

prepare_volcano_display_axes <- function(dt) {
  dt <- data.table::copy(dt)

  finite_abs_z <- abs(dt$z_plot[is.finite(dt$z_plot)])
  finite_logp <- dt$neglog10_p[is.finite(dt$neglog10_p)]

  z_cap <- if (length(finite_abs_z) > 0L) {
    min(12, max(8, as.numeric(stats::quantile(
      finite_abs_z,
      probs = 0.995,
      names = FALSE,
      na.rm = TRUE
    ))))
  } else {
    8
  }

  logp_cap <- if (length(finite_logp) > 0L) {
    min(30, max(10, as.numeric(stats::quantile(
      finite_logp,
      probs = 0.995,
      names = FALSE,
      na.rm = TRUE
    ))))
  } else {
    10
  }

  if (!is.finite(z_cap) || z_cap <= 0) {
    z_cap <- 8
  }
  if (!is.finite(logp_cap) || logp_cap <= 0) {
    logp_cap <- 10
  }

  dt[, z_display := pmax(pmin(z_plot, z_cap), -z_cap)]
  dt[, neglog10_p_display := pmin(neglog10_p, logp_cap)]

  list(
    table = dt,
    z_cap = z_cap,
    logp_cap = logp_cap
  )
}

create_score_plot <- function(
  bundle = NULL,
  table = NULL,
  which = NULL,
  test_tail = c("auto", "one_tailed", "two_tailed"),
  feature_column = NULL,
  z_column = "zstat",
  p_value_column = NULL,
  label_features = NULL,
  title = "conseguiR Scores"
) {
  required_packages <- c("ggplot2", "ggrepel", "ggnewscale")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, logical(1))]
  if (length(missing_packages) > 0L) {
    stop("The following packages are required to create score plots: ",
         paste(missing_packages, collapse = ", "))
  }

  if (is.null(table)) {
    if (is.null(bundle) || !is.list(bundle)) {
      stop("Provide either `bundle` or `table` to `create_score_plot()`.")
    }

    if (!is.null(which) && which %in% names(bundle)) {
      table <- bundle[[which]]
    } else if (!is.null(which) && "objects" %in% names(bundle) && which %in% names(bundle$objects)) {
      table <- bundle$objects[[which]]
    } else if ("objects" %in% names(bundle)) {
      candidates <- intersect(c("gene_scores", "reg_scores", "all_genes", "top_genes"), names(bundle$objects))
      if (length(candidates) == 0L) {
        stop("Could not resolve a plottable table from the supplied bundle.")
      }
      table <- bundle$objects[[candidates[[1]]]]
    } else {
      stop("Could not resolve a plottable table from the supplied bundle.")
    }
  }

  dt <- prepare_score_plot_table(
    table = table,
    which = which,
    feature_column = feature_column,
    z_column = z_column,
    p_value_column = p_value_column,
    label_features = label_features
  )

  plot_mode <- infer_plot_mode(dt, test_tail = test_tail)

  base <- ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#111827"),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (identical(plot_mode, "volcano")) {
    volcano_display <- prepare_volcano_display_axes(dt)
    dt <- volcano_display$table
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = z_display, y = neglog10_p_display)) +
      ggplot2::geom_point(
        data = dt[highlighted == FALSE],
        colour = "#b3b3b3",
        alpha = 0.8,
        size = 2
      ) +
      ggplot2::geom_point(
        data = dt[highlighted == TRUE],
        colour = "#8a1538",
        alpha = 0.75,
        size = 2.4
      ) +
      ggrepel::geom_text_repel(
        data = dt[should_label == TRUE],
        ggplot2::aes(label = ifelse(matched_label != "", matched_label, feature_label)),
        colour = "#111111",
        size = 3.3,
        box.padding = 0.35,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      ggplot2::labs(
        title = title,
        x = "Z-score",
        y = expression(-log[10](p-value))
      ) +
      ggplot2::coord_cartesian(
        xlim = c(-volcano_display$z_cap, volcano_display$z_cap),
        ylim = c(0, volcano_display$logp_cap)
      )
  } else {
    data.table::setorder(dt, -z_plot, feature_label)
    dt[, rank_plot := seq_len(.N)]
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = rank_plot, y = z_plot)) +
      ggplot2::geom_point(
        data = dt[highlighted == FALSE],
        colour = "#b3b3b3",
        alpha = 0.8,
        size = 2
      ) +
      ggplot2::geom_point(
        data = dt[highlighted == TRUE],
        colour = "#8a1538",
        alpha = 0.75,
        size = 2.4
      ) +
      ggrepel::geom_text_repel(
        data = dt[should_label == TRUE],
        ggplot2::aes(label = ifelse(matched_label != "", matched_label, feature_label)),
        colour = "#111111",
        size = 3.3,
        box.padding = 0.35,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf
      ) +
      ggplot2::labs(
        title = title,
        x = "Rank",
        y = "Z-score"
      )
  }

  list(
    plot = p + base,
    plot_data = dt,
    plot_mode = plot_mode
  )
}

save_score_plot <- function(
  bundle = NULL,
  table = NULL,
  file_path,
  which = NULL,
  test_tail = "auto",
  feature_column = NULL,
  z_column = "zstat",
  p_value_column = NULL,
  label_features = NULL,
  title = "conseguiR Scores",
  width = 10,
  height = 7,
  dpi = 300
) {
  if (missing(file_path) || is.null(file_path) || !nzchar(file_path)) {
    stop("`file_path` must be provided to save the score plot.")
  }

  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  out <- create_score_plot(
    bundle = bundle,
    table = table,
    which = which,
    test_tail = test_tail,
    feature_column = feature_column,
    z_column = z_column,
    p_value_column = p_value_column,
    label_features = label_features,
    title = title
  )

  ggplot2::ggsave(
    filename = file_path,
    plot = out$plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )

  invisible(file_path)
}

normalize_locus_chromosome <- function(chromosome) {
  chr <- as.character(chromosome)
  if (length(chr) == 0L) {
    stop("`chromosome` must be a non-empty string.")
  }
  chr <- trimws(chr)
  if (any(!nzchar(chr))) {
    stop("`chromosome` must be a non-empty string.")
  }
  has_prefix <- grepl("^chr", chr, ignore.case = TRUE)
  out <- chr
  out[has_prefix] <- paste0("chr", sub("^chr", "", chr[has_prefix], ignore.case = TRUE))
  out[!has_prefix] <- paste0("chr", chr[!has_prefix])
  out
}

assign_interval_lanes <- function(dt, start_col, end_col) {
  if (nrow(dt) == 0L) {
    dt[, lane := integer()]
    return(dt)
  }

  out <- data.table::copy(dt)
  out[, row_idx_internal := seq_len(.N)]
  data.table::setorderv(out, cols = c(start_col, end_col))

  lane_ends <- numeric()
  lane_idx <- integer(nrow(out))

  for (i in seq_len(nrow(out))) {
    start_i <- as.numeric(out[[start_col]][[i]])
    end_i <- as.numeric(out[[end_col]][[i]])
    lane <- which(start_i > lane_ends)[1]
    if (is.na(lane)) {
      lane <- length(lane_ends) + 1L
      lane_ends <- c(lane_ends, end_i)
    } else {
      lane_ends[[lane]] <- end_i
    }
    lane_idx[[i]] <- lane - 1L
  }

  out[, lane := lane_idx]
  data.table::setorderv(out, "row_idx_internal")
  out[, row_idx_internal := NULL]
  out
}

rescale_track_fill <- function(x) {
  x <- as.numeric(x)
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0L) {
    return(rep(0.15, length(x)))
  }

  xmin <- min(finite_x, na.rm = TRUE)
  xmax <- max(finite_x, na.rm = TRUE)
  if (!is.finite(xmin) || !is.finite(xmax) || identical(xmin, xmax)) {
    return(rep(0.8, length(x)))
  }

  scaled <- (x - xmin) / (xmax - xmin)
  scaled[!is.finite(scaled)] <- 0.15
  pmax(0.15, pmin(1, scaled))
}

safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[!is.finite(out)] <- NA_real_
  out
}

resolve_locus_selected_genes <- function(selected_subgraph = NULL, nodes = NULL, nodes_path = NULL) {
  selected_nodes <- nodes
  if (is.null(selected_nodes) && !is.null(selected_subgraph)) {
    if ("objects" %in% names(selected_subgraph) && "nodes" %in% names(selected_subgraph$objects)) {
      selected_nodes <- selected_subgraph$objects$nodes
    } else if ("nodes" %in% names(selected_subgraph)) {
      selected_nodes <- selected_subgraph$nodes
    }
  }
  if (is.null(selected_nodes) && !is.null(nodes_path) && file.exists(nodes_path)) {
    selected_nodes <- data.table::fread(nodes_path, showProgress = FALSE)
  }
  if (is.null(selected_nodes)) {
    return(character())
  }

  selected_nodes <- data.table::as.data.table(selected_nodes)
  gene_col <- intersect(c("gene_name", "node_id"), names(selected_nodes))
  if (length(gene_col) == 0L) {
    return(character())
  }

  unique(as.character(stats::na.omit(selected_nodes[[gene_col[[1]]]])))
}

read_minimal_locus_gwas_table <- function(gwas_sumstats = NULL) {
  if (is.null(gwas_sumstats)) {
    return(NULL)
  }

  if (is.character(gwas_sumstats) && length(gwas_sumstats) == 1L) {
    cache_key <- paste0("gwas_minimal::", normalizePath(gwas_sumstats, winslash = "/", mustWork = FALSE))
    if (exists(cache_key, envir = .conseguiR_plot_cache, inherits = FALSE)) {
      return(data.table::copy(get(cache_key, envir = .conseguiR_plot_cache, inherits = FALSE)))
    }

    header_cols <- names(data.table::fread(gwas_sumstats, nrows = 0L, showProgress = FALSE))
    snp_id_col <- intersect(c("hm_rsid", "rsid", "rs_id", "hm_variant_id", "variant_id"), header_cols)
    chr_col <- intersect(c("hm_chrom", "chromosome", "chr", "chrom"), header_cols)
    pos_col <- intersect(c("hm_pos", "base_pair_location", "position", "pos", "bp"), header_cols)
    p_col <- intersect(c("p_value", "pval", "p"), header_cols)
    if (length(snp_id_col) == 0L || length(chr_col) == 0L || length(pos_col) == 0L || length(p_col) == 0L) {
      return(NULL)
    }
    keep_cols <- unique(c(snp_id_col[[1]], chr_col[[1]], pos_col[[1]], p_col[[1]]))
    gwas_dt <- data.table::as.data.table(data.table::fread(gwas_sumstats, select = keep_cols, showProgress = FALSE))
  } else {
    gwas_dt <- data.table::as.data.table(data.table::copy(gwas_sumstats))
    snp_id_col <- intersect(c("hm_rsid", "rsid", "rs_id", "hm_variant_id", "variant_id"), names(gwas_dt))
    chr_col <- intersect(c("hm_chrom", "chromosome", "chr", "chrom"), names(gwas_dt))
    pos_col <- intersect(c("hm_pos", "base_pair_location", "position", "pos", "bp"), names(gwas_dt))
    p_col <- intersect(c("p_value", "pval", "p"), names(gwas_dt))
    if (length(snp_id_col) == 0L || length(chr_col) == 0L || length(pos_col) == 0L || length(p_col) == 0L) {
      return(NULL)
    }
    cache_key <- NULL
  }

  out <- data.table::copy(gwas_dt)
  out[, chromosome := normalize_locus_chromosome(get(chr_col[[1]]))]
  out[, position := suppressWarnings(as.integer(get(pos_col[[1]])))]
  out[, p_value := suppressWarnings(as.numeric(get(p_col[[1]])))]
  out[, rsid := tolower(as.character(get(snp_id_col[[1]])))]
  out <- out[
    !is.na(chromosome) & nzchar(chromosome) &
      !is.na(position) &
      !is.na(p_value) & is.finite(p_value) & p_value > 0 &
      !is.na(rsid) & nzchar(rsid),
    .(rsid, chromosome, position, p_value)
  ]
  out <- unique(out, by = c("rsid", "chromosome", "position", "p_value"))

  if (!is.null(cache_key)) {
    assign(cache_key, out, envir = .conseguiR_plot_cache)
  }

  data.table::copy(out)
}

prepare_locus_gwas_label <- function(gwas_sumstats = NULL, chromosome, start, end, reg_nodes = NULL) {
  if (is.null(gwas_sumstats)) {
    return(NULL)
  }

  target_start <- as.integer(start)
  target_end <- as.integer(end)
  target_reg_id <- NA_character_

  if (!is.null(reg_nodes)) {
    reg_dt <- data.table::as.data.table(data.table::copy(reg_nodes))
    reg_dt <- reg_dt[
      !is.na(germline_score) &
        is.finite(germline_score) &
        !is.na(feature_start) &
        !is.na(feature_end)
    ]
    if (nrow(reg_dt) > 0L) {
      data.table::setorderv(reg_dt, c("germline_score", "feature_start", "feature_end"), c(-1L, 1L, 1L))
      target_start <- as.integer(reg_dt$feature_start[[1]])
      target_end <- as.integer(reg_dt$feature_end[[1]])
      target_reg_id <- as.character(reg_dt$feature_id[[1]])
    }
  }

  gwas_dt <- read_minimal_locus_gwas_table(gwas_sumstats)
  if (is.null(gwas_dt) || nrow(gwas_dt) == 0L) {
    return(NULL)
  }

  target_chr <- normalize_locus_chromosome(chromosome)
  gwas_dt <- gwas_dt[
    chromosome == target_chr &
      position >= target_start &
      position <= target_end
  ]

  if (nrow(gwas_dt) == 0L) {
    return(NULL)
  }

  data.table::setorderv(gwas_dt, c("p_value", "position"), c(1L, 1L))
  top_snp <- gwas_dt[1]
  data.table::data.table(
    snp_id = top_snp$rsid[[1]],
    position = top_snp$position[[1]],
    p_value = top_snp$p_value[[1]],
    reg_feature_id = target_reg_id,
    reg_start = target_start,
    reg_end = target_end
  )
}

read_locus_gwas_hits <- function(gwas_sumstats = NULL, chromosome, start, end) {
  if (is.null(gwas_sumstats)) {
    return(NULL)
  }

  gwas_dt <- read_minimal_locus_gwas_table(gwas_sumstats)
  if (is.null(gwas_dt) || nrow(gwas_dt) == 0L) {
    return(NULL)
  }

  target_chr <- normalize_locus_chromosome(chromosome)
  gwas_dt <- gwas_dt[
    chromosome == target_chr &
      position >= as.integer(start) &
      position <= as.integer(end),
    .(rsid, position, p_value)
  ]

  if (nrow(gwas_dt) == 0L) {
    return(NULL)
  }

  unique(gwas_dt, by = c("rsid", "position", "p_value"))
}

fetch_dbsnp_rsid_pmids <- function(rsids, max_pmids_per_rsid = 200L, verbose = FALSE) {
  rsids <- unique(tolower(as.character(rsids %||% character())))
  rsids <- rsids[grepl("^rs[0-9]+$", rsids)]
  if (length(rsids) == 0L) {
    return(NULL)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package `jsonlite` is required to query NCBI E-utilities.")
  }

  max_pmids_per_rsid <- suppressWarnings(as.integer(max_pmids_per_rsid[[1]]))
  if (is.na(max_pmids_per_rsid) || max_pmids_per_rsid <= 0L) {
    max_pmids_per_rsid <- 200L
  }

  if (!exists("dbsnp_pmid_cache", envir = .conseguiR_plot_cache, inherits = FALSE)) {
    .conseguiR_plot_cache$dbsnp_pmid_cache <- new.env(parent = emptyenv())
  }
  cache_env <- .conseguiR_plot_cache$dbsnp_pmid_cache

  cached_mask <- vapply(
    rsids,
    FUN.VALUE = logical(1),
    FUN = function(rsid) exists(rsid, envir = cache_env, inherits = FALSE)
  )
  cached_rsids <- rsids[cached_mask]
  query_rsids <- rsids[!cached_mask]

  if (isTRUE(verbose)) {
    message(
      "Querying dbSNP publication pages for PMID-backed SNP evidence across ",
      length(query_rsids),
      " uncached rsID(s)",
      if (length(cached_rsids) > 0L) paste0(" (", length(cached_rsids), " served from cache)") else "",
      "."
    )
  }

  out_list <- vector("list", length(query_rsids))
  pb <- NULL
  if (isTRUE(verbose) && length(query_rsids) > 0L) {
    pb <- utils::txtProgressBar(min = 0, max = length(query_rsids), style = 3)
    on.exit({
      if (!is.null(pb)) {
        close(pb)
      }
    }, add = TRUE)
  }
  for (i in seq_along(query_rsids)) {
    rsid <- query_rsids[[i]]
    snp_id <- sub("^rs", "", rsid)
    endpoint <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=snp&db=pubmed&id=",
      snp_id,
      "&retmode=json"
    )
    payload <- tryCatch(jsonlite::fromJSON(endpoint, simplifyVector = FALSE), error = function(e) NULL)
    linksets <- payload$linksets %||% NULL
    pmid_hits <- character()
    if (!is.null(linksets) && length(linksets) > 0L) {
      for (j in seq_along(linksets)) {
        linkdbs <- linksets[[j]]$linksetdbs %||% NULL
        if (is.null(linkdbs) || length(linkdbs) == 0L) {
          next
        }
        for (k in seq_along(linkdbs)) {
          links <- linkdbs[[k]]$links %||% NULL
          if (!is.null(links)) {
            pmid_hits <- c(pmid_hits, as.character(links))
          }
        }
      }
    }
    pmid_hits <- unique(gsub("\\D", "", pmid_hits))
    pmid_hits <- pmid_hits[!is.na(pmid_hits) & nzchar(pmid_hits)]
    assign(rsid, pmid_hits, envir = cache_env)
    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, i)
    }
    if (length(pmid_hits) == 0L) {
      next
    }

    if (length(pmid_hits) > max_pmids_per_rsid) {
      pmid_hits <- pmid_hits[seq_len(max_pmids_per_rsid)]
    }

    out_list[[i]] <- data.table::data.table(rsid = rsid, pmid = pmid_hits, source = "dbsnp")
  }

  if (!is.null(pb) && length(query_rsids) > 0L) {
    utils::setTxtProgressBar(pb, length(query_rsids))
  }

  cached_list <- vector("list", length(cached_rsids))
  for (i in seq_along(cached_rsids)) {
    rsid <- cached_rsids[[i]]
    pmid_hits <- get(rsid, envir = cache_env, inherits = FALSE)
    if (length(pmid_hits) == 0L) {
      next
    }
    if (length(pmid_hits) > max_pmids_per_rsid) {
      pmid_hits <- pmid_hits[seq_len(max_pmids_per_rsid)]
    }
    cached_list[[i]] <- data.table::data.table(rsid = rsid, pmid = pmid_hits, source = "dbsnp")
  }

  out <- data.table::rbindlist(c(cached_list, out_list), use.names = TRUE, fill = TRUE)
  if (nrow(out) == 0L) {
    if (isTRUE(verbose)) {
      message("Found citation-backed evidence for 0 of ", length(rsids), " queried rsID(s).")
    }
    return(NULL)
  }

  if (isTRUE(verbose)) {
    message("Found citation-backed evidence for ", data.table::uniqueN(out$rsid), " of ", length(rsids), " queried rsID(s).")
  }

  unique(out)
}

prepare_staggered_snp_label_positions <- function(label_dt, locus_start, locus_end, base_y, lane_step = 0.28) {
  if (is.null(label_dt) || nrow(label_dt) == 0L) {
    return(label_dt)
  }

  out <- data.table::copy(label_dt)
  locus_width <- max(as.numeric(locus_end) - as.numeric(locus_start), 1)
  label_half_width <- max(25000, ceiling(locus_width * 0.018))
  out[, label_start := pmax(as.integer(position - label_half_width), as.integer(locus_start))]
  out[, label_end := pmin(as.integer(position + label_half_width), as.integer(locus_end))]
  lanes <- assign_interval_lanes(out[, .(position, label_start, label_end)], "label_start", "label_end")
  out[, label_lane := lanes$lane]
  out[, label_y := base_y + lane_step * label_lane]
  out[, label_x := as.numeric(position)]
  left_threshold <- as.numeric(locus_start) + label_half_width
  right_threshold <- as.numeric(locus_end) - label_half_width
  out[label_x < left_threshold, label_x := left_threshold]
  out[label_x > right_threshold, label_x := right_threshold]
  out[, c("label_start", "label_end") := NULL]
  out
}

fetch_variant_literature_pmids <- function(rsids, max_pmids_per_rsid = 200L, verbose = FALSE) {
  dbsnp_dt <- fetch_dbsnp_rsid_pmids(
    rsids = rsids,
    max_pmids_per_rsid = max_pmids_per_rsid,
    verbose = verbose
  )
  if (is.null(dbsnp_dt) || nrow(dbsnp_dt) == 0L) {
    return(NULL)
  }

  unique(dbsnp_dt)
}

normalize_requested_label_count <- function(top_n) {
  if (is.null(top_n) || length(top_n) == 0L) {
    return(0L)
  }

  if (is.infinite(top_n[[1]])) {
    return(.Machine$integer.max %/% 2L)
  }

  top_n <- suppressWarnings(as.integer(top_n[[1]]))
  if (is.na(top_n) || top_n <= 0L) {
    return(0L)
  }

  top_n
}

prepare_fallback_reg_snp_labels <- function(gwas_dt = NULL, reg_nodes = NULL, top_n = 3L) {
  top_n <- normalize_requested_label_count(top_n)
  if (top_n <= 0L || is.null(gwas_dt) || is.null(reg_nodes)) {
    return(NULL)
  }

  reg_dt <- data.table::as.data.table(data.table::copy(reg_nodes))[
    !is.na(germline_score) &
      is.finite(germline_score) &
      !is.na(feature_start) &
      !is.na(feature_end),
    .(
      reg_feature_id = as.character(feature_id),
      reg_start = as.integer(feature_start),
      reg_end = as.integer(feature_end),
      reg_germline_score = safe_numeric(germline_score)
    )
  ]
  if (nrow(reg_dt) == 0L) {
    return(NULL)
  }

  data.table::setorderv(reg_dt, c("reg_germline_score", "reg_start", "reg_end"), c(-1L, 1L, 1L))
  reg_dt <- reg_dt[seq_len(min(.N, top_n))]

  gwas_work <- data.table::copy(gwas_dt)
  gwas_work[, snp_start := position]
  gwas_work[, snp_end := position]
  data.table::setkey(reg_dt, reg_start, reg_end)
  overlaps <- data.table::foverlaps(
    gwas_work[, .(rsid, position, p_value, snp_start, snp_end)],
    reg_dt,
    by.x = c("snp_start", "snp_end"),
    by.y = c("reg_start", "reg_end"),
    nomatch = 0L
  )
  if (nrow(overlaps) == 0L) {
    return(NULL)
  }

  data.table::setorderv(overlaps, c("reg_feature_id", "p_value", "position"), c(1L, 1L, 1L))
  fallback <- overlaps[, .SD[1], by = reg_feature_id]
  data.table::setorderv(fallback, c("reg_germline_score", "p_value"), c(-1L, 1L))
  fallback[seq_len(min(.N, top_n))]
}

prepare_locus_lit_snp_labels <- function(
  gwas_sumstats = NULL,
  rsid_pmid = NULL,
  reg_nodes = NULL,
  chromosome,
  start,
  end,
  top_n = 3L,
  pmid_query = NULL,
  pmid_page_size = 1000L,
  verbose = FALSE
) {
  top_n <- normalize_requested_label_count(top_n)
  if (top_n <= 0L || is.null(gwas_sumstats) || is.null(reg_nodes)) {
    return(list(labels = NULL, mapping = NULL, mode = "none", requested_n = top_n, available_n = 0L))
  }

  gwas_dt <- read_locus_gwas_hits(
    gwas_sumstats = gwas_sumstats,
    chromosome = chromosome,
    start = start,
    end = end
  )
  if (is.null(gwas_dt) || nrow(gwas_dt) == 0L) {
    return(list(labels = NULL, mapping = NULL, mode = "none", requested_n = top_n, available_n = 0L))
  }

  reg_dt <- data.table::as.data.table(data.table::copy(reg_nodes))[
    !is.na(feature_start) & !is.na(feature_end),
    .(
      reg_feature_id = as.character(feature_id),
      reg_start = as.integer(feature_start),
      reg_end = as.integer(feature_end),
      reg_germline_score = safe_numeric(germline_score)
    )
  ]
  if (nrow(reg_dt) == 0L) {
    return(list(labels = NULL, mapping = NULL, mode = "none", requested_n = top_n, available_n = 0L))
  }

  gwas_dt[, snp_start := position]
  gwas_dt[, snp_end := position]
  gwas_dt[, rsid := tolower(rsid)]
  data.table::setkey(reg_dt, reg_start, reg_end)
  overlaps <- data.table::foverlaps(
    gwas_dt[, .(rsid, position, p_value, snp_start, snp_end)],
    reg_dt,
    by.x = c("snp_start", "snp_end"),
    by.y = c("reg_start", "reg_end"),
    nomatch = 0L
  )
  if (nrow(overlaps) == 0L) {
    fallback <- prepare_fallback_reg_snp_labels(gwas_dt = gwas_dt, reg_nodes = reg_nodes, top_n = top_n)
    return(list(
      labels = fallback,
      mapping = NULL,
      mode = if (is.null(fallback) || nrow(fallback) == 0L) "none" else "fallback",
      requested_n = top_n,
      available_n = if (is.null(fallback)) 0L else nrow(fallback)
    ))
  }

  all_overlap_rsids <- sort(unique(overlaps$rsid))
  if (isTRUE(verbose)) {
    message(
      "Querying citation support across ",
      length(all_overlap_rsids),
      " rsID(s) overlapping regulatory elements in this locus."
    )
  }

  pmid_dt <- NULL
  if (!is.null(rsid_pmid)) {
    if (is.character(rsid_pmid) && length(rsid_pmid) == 1L) {
      pmid_dt <- data.table::as.data.table(data.table::fread(rsid_pmid, showProgress = FALSE))
    } else {
      pmid_dt <- data.table::as.data.table(data.table::copy(rsid_pmid))
    }
  } else {
    pmid_dt <- fetch_variant_literature_pmids(
      rsids = all_overlap_rsids,
      max_pmids_per_rsid = pmid_page_size,
      verbose = verbose
    )
  }

  if (is.null(pmid_dt) || nrow(pmid_dt) == 0L) {
    fallback <- prepare_fallback_reg_snp_labels(gwas_dt = gwas_dt, reg_nodes = reg_nodes, top_n = top_n)
    return(list(
      labels = fallback,
      mapping = NULL,
      mode = if (is.null(fallback) || nrow(fallback) == 0L) "none" else "fallback",
      requested_n = top_n,
      available_n = if (is.null(fallback)) 0L else nrow(fallback)
    ))
  }

  rsid_col <- intersect(c("rsid", "hm_rsid", "rs_id"), names(pmid_dt))
  pmid_col <- intersect(c("pmid", "PMID"), names(pmid_dt))
  if (length(rsid_col) == 0L || length(pmid_col) == 0L) {
    fallback <- prepare_fallback_reg_snp_labels(gwas_dt = gwas_dt, reg_nodes = reg_nodes, top_n = top_n)
    return(list(
      labels = fallback,
      mapping = NULL,
      mode = if (is.null(fallback) || nrow(fallback) == 0L) "none" else "fallback",
      requested_n = top_n,
      available_n = if (is.null(fallback)) 0L else nrow(fallback)
    ))
  }

  pmid_dt <- pmid_dt[, .(
    rsid = tolower(as.character(get(rsid_col[[1]]))),
    pmid = as.character(get(pmid_col[[1]]))
  )]
  pmid_dt <- pmid_dt[!is.na(rsid) & nzchar(rsid) & !is.na(pmid) & nzchar(pmid)]
  if (nrow(pmid_dt) == 0L) {
    fallback <- prepare_fallback_reg_snp_labels(gwas_dt = gwas_dt, reg_nodes = reg_nodes, top_n = top_n)
    return(list(
      labels = fallback,
      mapping = NULL,
      mode = if (is.null(fallback) || nrow(fallback) == 0L) "none" else "fallback",
      requested_n = top_n,
      available_n = if (is.null(fallback)) 0L else nrow(fallback)
    ))
  }

  overlaps <- merge(overlaps, pmid_dt, by = "rsid", all = FALSE)
  if (nrow(overlaps) == 0L) {
    fallback <- prepare_fallback_reg_snp_labels(gwas_dt = gwas_dt, reg_nodes = reg_nodes, top_n = top_n)
    return(list(
      labels = fallback,
      mapping = NULL,
      mode = if (is.null(fallback) || nrow(fallback) == 0L) "none" else "fallback",
      requested_n = top_n,
      available_n = if (is.null(fallback)) 0L else nrow(fallback)
    ))
  }

  data.table::setorderv(overlaps, c("rsid", "reg_germline_score", "p_value"), c(1L, -1L, 1L))
  top_reg_per_snp <- unique(overlaps[, .(
    rsid,
    position,
    p_value,
    reg_feature_id,
    reg_start,
    reg_end,
    reg_germline_score
  )], by = "rsid")

  pmid_summary <- overlaps[, .(
    n_pmids = data.table::uniqueN(pmid),
    pmids = paste(sort(unique(pmid)), collapse = "|")
  ), by = rsid]

  label_dt <- merge(top_reg_per_snp, pmid_summary, by = "rsid", all.x = TRUE)
  data.table::setorderv(label_dt, c("p_value", "n_pmids", "reg_germline_score", "position"), c(1L, -1L, -1L, 1L))
  available_lit_n <- nrow(label_dt)
  label_dt <- label_dt[seq_len(min(.N, top_n))]
  if (nrow(label_dt) == 0L) {
    fallback <- prepare_fallback_reg_snp_labels(gwas_dt = gwas_dt, reg_nodes = reg_nodes, top_n = top_n)
    return(list(
      labels = fallback,
      mapping = NULL,
      mode = if (is.null(fallback) || nrow(fallback) == 0L) "none" else "fallback",
      requested_n = top_n,
      available_n = if (is.null(fallback)) 0L else nrow(fallback)
    ))
  }

  list(
    labels = data.table::copy(label_dt),
    mapping = data.table::copy(label_dt[, .(rsid, pmids, n_pmids, p_value, position, reg_feature_id, reg_start, reg_end, reg_germline_score)]),
    mode = "literature",
    requested_n = top_n,
    available_n = available_lit_n
  )
}

prepare_locus_plot_bundle <- function(
  nodes,
  edges,
  selected_genes = character(),
  chromosome,
  start,
  end,
  label_features = NULL,
  gwas_sumstats = NULL,
  label_top_gwas_snp = FALSE,
  rsid_pmid = NULL,
  label_top_lit_snps = 0L,
  pmid_query = NULL,
  pmid_page_size = 1000L,
  verbose = FALSE
) {
  locus_chr <- normalize_locus_chromosome(chromosome)
  locus_start <- as.integer(start)
  locus_end <- as.integer(end)
  if (is.na(locus_start) || is.na(locus_end) || locus_end < locus_start) {
    stop("`start` and `end` must define a valid genomic interval.")
  }

  nodes_dt <- data.table::as.data.table(data.table::copy(nodes))
  edges_dt <- data.table::as.data.table(data.table::copy(edges))
  gene_nodes <- nodes_dt[
    node_type == "gene" &
      chr == locus_chr &
      !is.na(start) & !is.na(end) &
      end >= locus_start & start <= locus_end
  ][, .(
    feature_id = as.character(node_id),
    feature_label = as.character(name),
    chromosome = as.character(chr),
    feature_start = as.integer(start),
    feature_end = as.integer(end),
    germline_score = safe_numeric(germline_score),
    somatic_score = safe_numeric(somatic_score),
    epigenomic_score = safe_numeric(epigenomic_score),
    post_germline = if ("post_germline" %in% names(.SD)) safe_numeric(post_germline) else NA_real_,
    post_somatic = if ("post_somatic" %in% names(.SD)) safe_numeric(post_somatic) else NA_real_,
    post_epigenomic = if ("post_epigenomic" %in% names(.SD)) safe_numeric(post_epigenomic) else NA_real_,
    post_norm = if ("post_norm" %in% names(.SD)) safe_numeric(post_norm) else NA_real_
  )]

  reg_nodes <- nodes_dt[
    node_type == "reg" &
      reg_chr == locus_chr &
      !is.na(reg_start) & !is.na(reg_end) &
      reg_end >= locus_start & reg_start <= locus_end
  ][, .(
    feature_id = as.character(node_id),
    feature_label = as.character(node_id),
    chromosome = as.character(reg_chr),
    feature_start = as.integer(reg_start),
    feature_end = as.integer(reg_end),
    germline_score = safe_numeric(germline_score),
    somatic_score = safe_numeric(somatic_score),
    epigenomic_score = safe_numeric(epigenomic_score)
  )]

  graph_gene_scores <- safe_numeric(nodes_dt[node_type == "gene", post_norm])
  gene_score_limits <- range(graph_gene_scores[is.finite(graph_gene_scores)], na.rm = TRUE)
  if (!all(is.finite(gene_score_limits))) {
    gene_score_limits <- c(0, 1)
  } else if (diff(gene_score_limits) <= 0) {
    gene_score_limits <- c(gene_score_limits[[1]], gene_score_limits[[1]] + 1)
  }

  gwas_label_dt <- if (isTRUE(label_top_gwas_snp)) {
    prepare_locus_gwas_label(
      gwas_sumstats = gwas_sumstats,
      chromosome = locus_chr,
      start = locus_start,
      end = locus_end,
      reg_nodes = reg_nodes
    )
  } else {
    NULL
  }

  lit_snp_info <- prepare_locus_lit_snp_labels(
    gwas_sumstats = gwas_sumstats,
    rsid_pmid = rsid_pmid,
    reg_nodes = reg_nodes,
    chromosome = locus_chr,
    start = locus_start,
    end = locus_end,
    top_n = label_top_lit_snps,
    pmid_query = pmid_query,
    pmid_page_size = pmid_page_size,
    verbose = verbose
  )

  gwas_hits_dt <- read_locus_gwas_hits(
    gwas_sumstats = gwas_sumstats,
    chromosome = locus_chr,
    start = locus_start,
    end = locus_end
  )
  gwas_logp_limits <- c(0, 1)
  if (!is.null(gwas_hits_dt) && nrow(gwas_hits_dt) > 0L) {
    gwas_hits_dt[, logp := -log10(pmax(p_value, 1e-300))]
    finite_logp <- gwas_hits_dt$logp[is.finite(gwas_hits_dt$logp)]
    if (length(finite_logp) > 0L) {
      gwas_logp_limits <- c(0, max(finite_logp))
      if (diff(gwas_logp_limits) <= 0) {
        gwas_logp_limits[[2]] <- gwas_logp_limits[[1]] + 1
      }
    }
  }

  reg_nodes[, current_norm := sqrt(
    fifelse(is.na(germline_score), 0, germline_score)^2 +
      fifelse(is.na(somatic_score), 0, somatic_score)^2 +
      fifelse(is.na(epigenomic_score), 0, epigenomic_score)^2
  )]

  selected_gene_set <- unique(as.character(c(selected_genes, label_features %||% character())))
  selected_gene_set <- selected_gene_set[nzchar(selected_gene_set)]
  gene_nodes[, highlighted := feature_label %in% selected_gene_set]

  reg_label_map <- read_backend_reg_label_map()
  if (nrow(reg_label_map) > 0L && nrow(reg_nodes) > 0L) {
    reg_nodes <- merge(reg_nodes, reg_label_map, by.x = "feature_id", by.y = "reg_elem_id", all.x = TRUE)
    reg_nodes[, linked_label := fifelse(!is.na(label) & label != "", label, feature_label)]
    reg_nodes[, label := NULL]
  } else {
    reg_nodes[, linked_label := feature_label]
  }

  track_defs <- data.table::data.table(
    track_name = c(
      "Reg element somatic z",
      "Reg element epigenomic z",
      "Reg element germline z",
      "GWAS locus SNPs",
      "Reg elements",
      "Genes (post-diffusion)"
    ),
    track_id = 6:1
  )

  gene_long <- data.table::rbindlist(list(
    gene_nodes[, .(feature_type = "gene", feature_id, feature_label, feature_start, feature_end, track_name = "Genes (post-diffusion)", score = post_norm, linked_label = feature_label, highlighted)]
  ), use.names = TRUE)

  reg_long <- data.table::rbindlist(list(
    reg_nodes[, .(feature_type = "reg", feature_id, feature_label, feature_start, feature_end, track_name = "Reg element somatic z", score = somatic_score, linked_label, highlighted = FALSE)],
    reg_nodes[, .(feature_type = "reg", feature_id, feature_label, feature_start, feature_end, track_name = "Reg element epigenomic z", score = epigenomic_score, linked_label, highlighted = FALSE)],
    reg_nodes[, .(feature_type = "reg", feature_id, feature_label, feature_start, feature_end, track_name = "Reg element germline z", score = germline_score, linked_label, highlighted = FALSE)],
    reg_nodes[, .(feature_type = "reg", feature_id, feature_label, feature_start, feature_end, track_name = "Reg elements", score = current_norm, linked_label, highlighted = FALSE)]
  ), use.names = TRUE)

  gwas_long <- if (!is.null(gwas_hits_dt) && nrow(gwas_hits_dt) > 0L) {
    gwas_hits_dt[, .(
      feature_type = "snp",
      feature_id = rsid,
      feature_label = rsid,
      feature_start = as.integer(position),
      feature_end = as.integer(position),
      track_name = "GWAS locus SNPs",
      score = logp,
      linked_label = rsid,
      highlighted = FALSE
    )]
  } else {
    NULL
  }

  plot_dt <- data.table::rbindlist(list(gene_long, reg_long, gwas_long), use.names = TRUE, fill = TRUE)
  plot_dt <- merge(plot_dt, track_defs, by = "track_name", all.x = TRUE)
  plot_dt[, feature_mid := (feature_start + feature_end) / 2]
  plot_dt <- plot_dt[feature_end >= locus_start & feature_start <= locus_end]

  gene_lanes <- assign_interval_lanes(unique(gene_nodes[, .(feature_id, feature_start, feature_end)]), "feature_start", "feature_end")
  reg_lanes <- assign_interval_lanes(unique(reg_nodes[, .(feature_id, feature_start, feature_end)]), "feature_start", "feature_end")
  plot_dt <- merge(plot_dt, gene_lanes[, .(feature_id, gene_lane = lane)], by = "feature_id", all.x = TRUE)
  plot_dt <- merge(plot_dt, reg_lanes[, .(feature_id, reg_lane = lane)], by = "feature_id", all.x = TRUE)
  plot_dt[is.na(gene_lane), gene_lane := 0L]
  plot_dt[is.na(reg_lane), reg_lane := 0L]
  plot_dt[, y := fifelse(
    feature_type == "gene",
    track_id - 0.05 - 0.18 * gene_lane,
    fifelse(
      feature_type == "reg",
      track_id - 0.02 - 0.18 * reg_lane,
      track_id + 0.02 + 0.56 * (score / max(gwas_logp_limits[[2]], 1))
    )
  )]

  reg_score_limits <- range(plot_dt[feature_type == "reg" & track_name != "Reg elements", score], na.rm = TRUE)
  if (!all(is.finite(reg_score_limits))) {
    reg_score_limits <- c(-1, 1)
  } else {
    reg_abs_max <- max(abs(reg_score_limits))
    reg_score_limits <- c(-reg_abs_max, reg_abs_max)
  }

  reg_norm_limits <- range(reg_nodes$current_norm[is.finite(reg_nodes$current_norm)], na.rm = TRUE)
  if (!all(is.finite(reg_norm_limits))) {
    reg_norm_limits <- c(0, 1)
  } else if (diff(reg_norm_limits) <= 0) {
    reg_norm_limits <- c(reg_norm_limits[[1]], reg_norm_limits[[1]] + 1)
  }

  edges_locus <- edges_dt[
    reg_chr == locus_chr & gene_chr == locus_chr &
      reg_end >= locus_start & reg_start <= locus_end &
      gene_end >= locus_start & gene_start <= locus_end
  ][, .(
    gene_id = as.character(from),
    reg_id = as.character(to),
    reg_start = as.integer(reg_start),
    reg_end = as.integer(reg_end),
    gene_start = as.integer(gene_start),
    gene_end = as.integer(gene_end),
    link_score = safe_numeric(link_score),
    confidence = safe_numeric(confidence),
    link_method = as.character(link_method)
  )]

  gene_lookup <- unique(gene_nodes[, .(gene_id = feature_id, gene_label = feature_label)])
  reg_lookup <- unique(reg_nodes[, .(reg_id = feature_id, linked_label, current_norm)])
  edges_locus <- merge(edges_locus, gene_lookup, by = "gene_id", all.x = TRUE)
  edges_locus <- merge(edges_locus, reg_lookup, by = "reg_id", all.x = TRUE)
  edges_locus[, gene_mid := (gene_start + gene_end) / 2]
  edges_locus[, reg_mid := (reg_start + reg_end) / 2]
  edges_locus[, interaction_strength := fifelse(
    is.finite(confidence),
    confidence,
    fifelse(is.finite(link_score), link_score, NA_real_)
  )]
  finite_strength <- edges_locus$interaction_strength[is.finite(edges_locus$interaction_strength)]
  if (length(finite_strength) > 0L) {
    strength_min <- min(finite_strength, na.rm = TRUE)
    strength_cap <- as.numeric(stats::quantile(finite_strength, probs = 0.95, na.rm = TRUE, names = FALSE, type = 7))
    if (!is.finite(strength_cap)) {
      strength_cap <- max(finite_strength, na.rm = TRUE)
    }
    strength_max <- max(strength_min, strength_cap)
    if (isTRUE(all.equal(strength_min, strength_max))) {
      edges_locus[, link_alpha := 0.65]
    } else {
      edges_locus[, strength_scaled := log1p(pmin(interaction_strength, strength_max))]
      scaled_min <- min(edges_locus$strength_scaled[is.finite(edges_locus$strength_scaled)], na.rm = TRUE)
      scaled_max <- max(edges_locus$strength_scaled[is.finite(edges_locus$strength_scaled)], na.rm = TRUE)
      edges_locus[, link_alpha := 0.15 + 0.75 * ((strength_scaled - scaled_min) / max(scaled_max - scaled_min, 1e-8))]
      edges_locus[!is.finite(link_alpha), link_alpha := 0.15]
      edges_locus[, strength_scaled := NULL]
    }
  } else {
    edges_locus[, link_alpha := 0.25]
  }

  if (length(selected_gene_set) > 0L) {
    edges_locus <- edges_locus[gene_label %in% selected_gene_set | vapply(
      strsplit(linked_label %||% "", "|", fixed = TRUE),
      function(x) any(x %in% selected_gene_set),
      logical(1)
    )]
  }

  reg_label_dt <- unique(reg_nodes[, .(feature_id, linked_label, current_norm, feature_start, feature_end)])
  reg_label_dt[, matched_label := vapply(
    strsplit(linked_label, "|", fixed = TRUE),
    function(x) {
      hit <- intersect(x, selected_gene_set)
      if (length(hit) == 0L) "" else hit[[1]]
    },
    character(1)
  )]
  reg_label_dt <- reg_label_dt[matched_label != ""]
  if (nrow(reg_label_dt) > 0L) {
    data.table::setorderv(reg_label_dt, c("matched_label", "current_norm"), order = c(1L, -1L))
    reg_label_dt <- reg_label_dt[, .SD[1], by = matched_label]
    reg_label_dt[, reg_mid := (feature_start + feature_end) / 2]
  }

  list(
    tracks = track_defs,
    features = plot_dt,
    gene_nodes = gene_nodes,
    reg_nodes = reg_nodes,
    links = edges_locus,
    reg_labels = reg_label_dt,
    gwas_label = gwas_label_dt,
    lit_snp_labels = lit_snp_info$labels,
    lit_snp_mapping = lit_snp_info$mapping,
    snp_label_mode = lit_snp_info$mode,
    snp_label_requested_n = lit_snp_info$requested_n %||% 0L,
    snp_label_available_n = lit_snp_info$available_n %||% 0L,
    snp_label_labeled_n = if (is.null(lit_snp_info$labels)) 0L else nrow(lit_snp_info$labels),
    selected_gene_set = selected_gene_set,
    locus = list(chromosome = locus_chr, start = locus_start, end = locus_end),
    gene_score_limits = gene_score_limits,
    reg_score_limits = reg_score_limits,
    reg_norm_limits = reg_norm_limits,
    gwas_logp_limits = gwas_logp_limits
  )
}

create_locus_context_plot <- function(
  nodes,
  edges,
  selected_genes = character(),
  chromosome,
  start,
  end,
  label_features = NULL,
  title = NULL,
  gwas_sumstats = NULL,
  label_top_gwas_snp = FALSE,
  rsid_pmid = NULL,
  label_top_lit_snps = 0L,
  pmid_query = NULL,
  pmid_page_size = 1000L,
  verbose = FALSE
) {
  required_packages <- c("ggplot2", "ggnewscale")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, logical(1))]
  if (length(missing_packages) > 0L) {
    stop("The following packages are required to create locus plots: ",
         paste(missing_packages, collapse = ", "))
  }

  bundle <- prepare_locus_plot_bundle(
    nodes = nodes,
    edges = edges,
    selected_genes = selected_genes,
    chromosome = chromosome,
    start = start,
    end = end,
    label_features = label_features,
    gwas_sumstats = gwas_sumstats,
    label_top_gwas_snp = label_top_gwas_snp,
    rsid_pmid = rsid_pmid,
    label_top_lit_snps = label_top_lit_snps,
    pmid_query = pmid_query,
    pmid_page_size = pmid_page_size,
    verbose = verbose
  )

  locus <- bundle$locus
  title <- title %||% paste0("Locus context: ", locus$chromosome, ":", format(locus$start, big.mark = ","), "-", format(locus$end, big.mark = ","))

  features_dt <- bundle$features
  links_dt <- bundle$links
  gene_label_dt <- unique(features_dt[feature_type == "gene", .(
    feature_label,
    feature_mid,
    y
  )])
  if (nrow(gene_label_dt) > 0L) {
    gene_label_dt[, y_label := y - 0.30]
  }

  gene_width <- max(ceiling((locus$end - locus$start) * 0.025), 60000L)
  features_dt[feature_type == "gene", `:=`(
    feature_start = as.integer(round(feature_mid - gene_width / 2)),
    feature_end = as.integer(round(feature_mid + gene_width / 2))
  )]
  features_dt[feature_type == "gene", xmid := (feature_start + feature_end) / 2]
  features_dt[feature_type == "reg", xmid := (feature_start + feature_end) / 2]
  features_dt[feature_type == "snp", xmid := feature_start]
  reg_track_id <- bundle$tracks[track_name == "Reg elements", track_id][[1]]
  gene_track_id <- bundle$tracks[track_name == "Genes (post-diffusion)", track_id][[1]]
  links_dt[, reg_y := reg_track_id - 0.02]
  links_dt[, gene_y := gene_track_id - 0.05]

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(
      data = bundle$tracks,
      ggplot2::aes(yintercept = track_id),
      colour = "#e5e7eb",
      linewidth = 0.4
    )

  if (nrow(features_dt[feature_type == "snp"]) > 0L) {
    p <- p + ggplot2::geom_point(
      data = features_dt[feature_type == "snp"],
      ggplot2::aes(x = xmid, y = y),
      inherit.aes = FALSE,
      shape = 16,
      colour = "#334155",
      alpha = 0.85,
      size = 1.7,
      show.legend = FALSE
    )
  }

  if (nrow(links_dt) > 0L) {
    p <- p + ggplot2::geom_curve(
      data = links_dt,
      ggplot2::aes(
        x = reg_mid,
        y = reg_y,
        xend = gene_mid,
        yend = gene_y
      ),
      curvature = -0.18,
      colour = "#111111",
      alpha = links_dt$link_alpha,
      linewidth = 0.25,
      show.legend = FALSE,
      inherit.aes = FALSE
    )
  }

  p <- p +
    ggplot2::geom_point(
      data = features_dt[feature_type == "reg" & track_name != "Reg elements"],
      ggplot2::aes(
        x = xmid,
        y = y,
        fill = score
      ),
      shape = 21,
      colour = "#4b5563",
      size = 3.6,
      stroke = 0.25,
      show.legend = TRUE
    ) +
    ggplot2::geom_point(
      data = features_dt[feature_type == "reg" & track_name == "Reg elements"],
      ggplot2::aes(
        x = xmid,
        y = y
      ),
      shape = 21,
      fill = "white",
      colour = "#4b5563",
      size = 3.6,
      stroke = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2563eb",
      mid = "#f8fafc",
      high = "#b91c1c",
      midpoint = 0,
      limits = bundle$reg_score_limits,
      name = "Reg element\nz-score",
      oob = scales::squish,
      guide = ggplot2::guide_colorbar(order = 1)
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(
      data = features_dt[feature_type == "reg" & track_name == "Reg elements"],
      ggplot2::aes(
        x = xmid,
        y = y,
        fill = score
      ),
      shape = 21,
      colour = "#4b5563",
      size = 3.6,
      stroke = 0.4,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#fee2e2", "#ef4444", "#b91c1c", "#7f1d1d"),
      limits = bundle$reg_norm_limits,
      name = "Reg element\ncombined score",
      oob = scales::squish,
      guide = ggplot2::guide_colorbar(order = 2)
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(
      data = features_dt[feature_type == "gene"],
      ggplot2::aes(
        xmin = feature_start,
        xmax = feature_end,
        ymin = y - 0.15,
        ymax = y + 0.15,
        fill = score
      ),
      colour = "#7f1d1d",
      linewidth = 0.55,
      show.legend = TRUE
    )

  if (nrow(gene_label_dt) > 0L) {
    p <- p + ggplot2::geom_text(
      data = gene_label_dt,
      ggplot2::aes(x = feature_mid, y = y_label, label = feature_label),
      inherit.aes = FALSE,
      colour = "#111111",
      size = 3.0,
      vjust = 1
    )
  }

  if (!is.null(bundle$lit_snp_labels) && nrow(bundle$lit_snp_labels) > 0L) {
    snp_y <- max(bundle$tracks$track_id, na.rm = TRUE) + 0.42
    lit_label_dt <- prepare_staggered_snp_label_positions(
      label_dt = bundle$lit_snp_labels,
      locus_start = locus$start,
      locus_end = locus$end,
      base_y = snp_y
    )
    p <- p +
      ggplot2::geom_vline(
        data = lit_label_dt,
        ggplot2::aes(xintercept = position),
        colour = "#0f172a",
        linewidth = 0.35,
        linetype = "22",
        alpha = 0.8,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      ggplot2::geom_label(
        data = lit_label_dt,
        ggplot2::aes(x = label_x, y = label_y, label = rsid),
        inherit.aes = FALSE,
        size = 2.8,
        linewidth = 0.2,
        label.padding = grid::unit(0.12, "lines"),
        fill = "#ffffff",
        colour = "#111111"
      )
  } else if (!is.null(bundle$gwas_label) && nrow(bundle$gwas_label) > 0L) {
    snp_y <- max(bundle$tracks$track_id, na.rm = TRUE) + 0.42
    gwas_label_dt <- prepare_staggered_snp_label_positions(
      label_dt = bundle$gwas_label[, .(position, rsid = snp_id)],
      locus_start = locus$start,
      locus_end = locus$end,
      base_y = snp_y
    )
    p <- p +
      ggplot2::geom_vline(
        data = gwas_label_dt,
        ggplot2::aes(xintercept = position),
        colour = "#0f172a",
        linewidth = 0.35,
        linetype = "22",
        alpha = 0.8,
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      ggplot2::geom_label(
        data = gwas_label_dt,
        ggplot2::aes(x = label_x, y = label_y, label = rsid),
        inherit.aes = FALSE,
        size = 2.8,
        linewidth = 0.2,
        label.padding = grid::unit(0.12, "lines"),
        fill = "#ffffff",
        colour = "#111111"
      )
  }

  snp_label_top_y <- NULL
  if (exists("lit_label_dt", inherits = FALSE) && !is.null(lit_label_dt) && nrow(lit_label_dt) > 0L) {
    snp_label_top_y <- max(lit_label_dt$label_y, na.rm = TRUE)
  } else if (exists("gwas_label_dt", inherits = FALSE) && !is.null(gwas_label_dt) && nrow(gwas_label_dt) > 0L) {
    snp_label_top_y <- max(gwas_label_dt$label_y, na.rm = TRUE)
  }

  top_ylim <- max(bundle$tracks$track_id) + 0.8
  if (!is.null(snp_label_top_y) && is.finite(snp_label_top_y)) {
    top_ylim <- max(top_ylim, snp_label_top_y + 0.35)
  }

  p <- p +
    ggplot2::scale_fill_gradientn(
      colours = c("#fee2e2", "#ef4444", "#b91c1c", "#7f1d1d"),
      limits = bundle$gene_score_limits,
      name = "Post-diffusion\nconseguiR score",
      oob = scales::squish,
      guide = ggplot2::guide_colorbar(order = 3)
    ) +
    ggplot2::scale_y_continuous(
      breaks = bundle$tracks$track_id,
      labels = bundle$tracks$track_name,
      minor_breaks = NULL
    ) +
    ggplot2::coord_cartesian(
      xlim = c(locus$start, locus$end),
      ylim = c(min(bundle$tracks$track_id) - 1.4, top_ylim),
      expand = FALSE,
      clip = "off"
    ) +
    ggplot2::labs(
      title = title,
      x = paste0(locus$chromosome, " position"),
      y = NULL
    ) +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", colour = "#111827"),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(face = "bold", colour = "#111827"),
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = grid::unit(c(10, 18, 10, 10), "pt")
    )

  list(
    plot = p,
    plot_data = bundle
  )
}

save_locus_context_plot <- function(
  nodes,
  edges,
  selected_genes = character(),
  chromosome,
  start,
  end,
  file_path,
  label_features = NULL,
  title = NULL,
  gwas_sumstats = NULL,
  label_top_gwas_snp = FALSE,
  rsid_pmid = NULL,
  label_top_lit_snps = 0L,
  pmid_query = NULL,
  pmid_page_size = 1000L,
  width = 14,
  height = 9,
  dpi = 300
) {
  if (missing(file_path) || is.null(file_path) || !nzchar(file_path)) {
    stop("`file_path` must be provided to save the locus plot.")
  }

  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  out <- create_locus_context_plot(
    nodes = nodes,
    edges = edges,
    selected_genes = selected_genes,
    chromosome = chromosome,
    start = start,
    end = end,
    label_features = label_features,
    title = title,
    gwas_sumstats = gwas_sumstats,
    label_top_gwas_snp = label_top_gwas_snp,
    rsid_pmid = rsid_pmid,
    label_top_lit_snps = label_top_lit_snps,
    pmid_query = pmid_query,
    pmid_page_size = pmid_page_size
  )

  ggplot2::ggsave(
    filename = file_path,
    plot = out$plot,
    width = width,
    height = height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )

  invisible(file_path)
}
