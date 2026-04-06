#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

conseguiR_plot_runtime_file <- function(relpath) {
  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(candidate)
  }

  pkg_path <- system.file(relpath, package = "conseguiR")
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }

  stop("Could not locate required plotting runtime file: ", relpath)
}

default_stage_plot_config <- list(
  score_output_prefix = "data/processed/conseguiR_score_plot"
)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

read_backend_gene_label_map <- function() {
  loc_path <- conseguiR_plot_runtime_file("inst/extdata/backend/NCBI38.gene.loc")
  dt <- data.table::fread(loc_path, header = FALSE, showProgress = FALSE)
  if (ncol(dt) < 6L) {
    stop("Backend gene location file does not contain the expected label column.")
  }

  out <- unique(dt[, .(
    feature_id = as.character(V1),
    label = as.character(V6)
  )], by = "feature_id")
  out[!is.na(feature_id) & feature_id != "" & !is.na(label) & label != ""]
}

read_backend_reg_label_map <- function() {
  mapping_path <- conseguiR_plot_runtime_file("inst/extdata/backend/genehancer_reg_target_labels.tsv.gz")
  if (!file.exists(mapping_path)) {
    return(data.table(reg_elem_id = character(), label = character()))
  }

  dt <- data.table::as.data.table(data.table::fread(mapping_path, showProgress = FALSE))
  unique(dt[!is.na(reg_elem_id) & reg_elem_id != "" & !is.na(label) & label != ""], by = "reg_elem_id")
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
  required_packages <- c("ggplot2", "ggrepel")
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
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = z_plot, y = neglog10_p)) +
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
