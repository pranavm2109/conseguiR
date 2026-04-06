#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

default_stage_plot_config <- list(
  scores_output_prefix = "data/processed/conseguiR_scores_plot",
  diffusion_output_prefix = "data/processed/conseguiR_diffusion_plot"
)

infer_feature_column <- function(dt, candidates = NULL) {
  dt <- data.table::as.data.table(dt)
  candidates <- candidates %||% c(
    "gene_name", "feature_id", "node_id", "reg_id", "symbol", "gene", "feature"
  )
  matches <- intersect(candidates, names(dt))
  if (length(matches) > 0L) {
    return(matches[[1]])
  }

  char_cols <- names(dt)[vapply(dt, function(x) is.character(x) || is.factor(x), logical(1))]
  if (length(char_cols) > 0L) {
    return(char_cols[[1]])
  }

  stop("Could not infer a feature-label column for plotting.")
}

infer_value_column <- function(dt, candidates = NULL) {
  dt <- data.table::as.data.table(dt)
  candidates <- candidates %||% c(
    "post_norm", "zscore", "z", "ZSTAT", "prize", "score",
    "wmis_cv", "qglobal_cv", "n_muts", "mean_signal_z"
  )
  matches <- intersect(candidates, names(dt))
  if (length(matches) > 0L) {
    return(matches[[1]])
  }

  numeric_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  numeric_cols <- setdiff(
    numeric_cols,
    c("chr", "chromosome", "start", "end", "pos", "position", "n", "rank")
  )
  if (length(numeric_cols) > 0L) {
    return(numeric_cols[[1]])
  }

  stop("Could not infer a numeric score column for plotting.")
}

prepare_ranked_plot_table <- function(
  table,
  feature_column = NULL,
  value_column = NULL,
  top_n = NULL,
  highlight_features = NULL
) {
  dt <- data.table::as.data.table(data.table::copy(table))
  feature_column <- feature_column %||% infer_feature_column(dt)
  value_column <- value_column %||% infer_value_column(dt)

  dt <- dt[!is.na(get(feature_column)) & !is.na(get(value_column))]
  dt[, feature := as.character(get(feature_column))]
  dt[, value := as.numeric(get(value_column))]
  dt <- dt[is.finite(value)]

  if (nrow(dt) == 0L) {
    stop("No finite rows are available to plot.")
  }

  data.table::setorder(dt, -value, feature)
  dt[, rank := seq_len(.N)]
  dt[, highlighted := feature %in% (highlight_features %||% character())]

  if (!is.null(top_n)) {
    dt <- dt[seq_len(min(as.integer(top_n), .N))]
  }

  dt[]
}

create_scores_plot <- function(
  bundle = NULL,
  table = NULL,
  which = NULL,
  plot_type = c("ranked_points", "top_bar", "histogram"),
  top_n = 25L,
  feature_column = NULL,
  value_column = NULL,
  highlight_features = NULL,
  title = "conseguiR Scores",
  subtitle = NULL
) {
  required_packages <- c("ggplot2", "scales")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, logical(1))]
  if (length(missing_packages) > 0L) {
    stop("The following packages are required to create score plots: ",
         paste(missing_packages, collapse = ", "))
  }

  plot_type <- match.arg(plot_type)

  if (is.null(table)) {
    if (is.null(bundle) || !is.list(bundle)) {
      stop("Provide either `bundle` or `table` to `create_scores_plot()`.")
    }

    if (!is.null(which) && which %in% names(bundle)) {
      table <- bundle[[which]]
    } else if (!is.null(which) && "objects" %in% names(bundle) && which %in% names(bundle$objects)) {
      table <- bundle$objects[[which]]
    } else if ("objects" %in% names(bundle)) {
      score_names <- intersect(
        c("gene_scores", "reg_scores", "scores"),
        names(bundle$objects)
      )
      if (length(score_names) == 0L) {
        stop("Could not resolve a score table from the supplied bundle.")
      }
      table <- bundle$objects[[score_names[[1]]]]
    } else {
      stop("Could not resolve a score table from the supplied bundle.")
    }
  }

  ranked_dt <- prepare_ranked_plot_table(
    table = table,
    feature_column = feature_column,
    value_column = value_column,
    top_n = if (plot_type %in% c("ranked_points", "top_bar")) top_n else NULL,
    highlight_features = highlight_features
  )

  value_label <- value_column %||% infer_value_column(data.table::as.data.table(table))
  p <- switch(
    plot_type,
    ranked_points = {
      ggplot2::ggplot(ranked_dt, ggplot2::aes(x = rank, y = value)) +
        ggplot2::geom_line(colour = "#8a1538", linewidth = 0.45, alpha = 0.6) +
        ggplot2::geom_point(
          ggplot2::aes(colour = highlighted),
          size = 2.2,
          show.legend = FALSE
        ) +
        ggplot2::scale_colour_manual(values = c(`TRUE` = "#b91c1c", `FALSE` = "#1f4e79")) +
        ggplot2::labs(x = "Rank", y = value_label, title = title, subtitle = subtitle)
    },
    top_bar = {
      ranked_dt[, feature := factor(feature, levels = rev(feature))]
      ggplot2::ggplot(ranked_dt, ggplot2::aes(x = feature, y = value, fill = highlighted)) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c(`TRUE` = "#b91c1c", `FALSE` = "#8a1538")) +
        ggplot2::labs(x = NULL, y = value_label, title = title, subtitle = subtitle)
    },
    histogram = {
      full_dt <- prepare_ranked_plot_table(
        table = table,
        feature_column = feature_column,
        value_column = value_column,
        top_n = NULL,
        highlight_features = highlight_features
      )
      ggplot2::ggplot(full_dt, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(fill = "#8a1538", colour = "white", bins = 30) +
        ggplot2::labs(x = value_label, y = "Count", title = title, subtitle = subtitle)
    }
  )

  p +
    ggplot2::theme_minimal(base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#111827"),
      plot.subtitle = ggplot2::element_text(colour = "#4b5563"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}

save_scores_plot <- function(
  bundle = NULL,
  table = NULL,
  file_path,
  which = NULL,
  plot_type = "ranked_points",
  top_n = 25L,
  feature_column = NULL,
  value_column = NULL,
  highlight_features = NULL,
  title = "conseguiR Scores",
  subtitle = NULL,
  width = 10,
  height = 7,
  dpi = 300
) {
  if (missing(file_path) || is.null(file_path) || !nzchar(file_path)) {
    stop("`file_path` must be provided to save the score plot.")
  }

  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  plot_obj <- create_scores_plot(
    bundle = bundle,
    table = table,
    which = which,
    plot_type = plot_type,
    top_n = top_n,
    feature_column = feature_column,
    value_column = value_column,
    highlight_features = highlight_features,
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

create_diffusion_plot <- function(
  diffusion = NULL,
  table = NULL,
  which = c("all_genes", "top_genes"),
  plot_type = c("ranked_points", "top_bar", "histogram"),
  top_n = 50L,
  gene_column = NULL,
  score_column = NULL,
  highlight_genes = NULL,
  title = "conseguiR Diffusion Scores",
  subtitle = NULL
) {
  which <- match.arg(which)
  plot_type <- match.arg(plot_type)

  if (is.null(table)) {
    if (is.null(diffusion) || !is.list(diffusion)) {
      stop("Provide either `diffusion` or `table` to `create_diffusion_plot()`.")
    }
    if (which %in% names(diffusion)) {
      table <- diffusion[[which]]
    } else if ("objects" %in% names(diffusion) && which %in% names(diffusion$objects)) {
      table <- diffusion$objects[[which]]
    } else {
      stop("Could not resolve the requested diffusion table: ", which)
    }
  }

  create_scores_plot(
    table = table,
    plot_type = plot_type,
    top_n = top_n,
    feature_column = gene_column %||% "gene_name",
    value_column = score_column,
    highlight_features = highlight_genes,
    title = title,
    subtitle = subtitle
  )
}

save_diffusion_plot <- function(
  diffusion = NULL,
  table = NULL,
  file_path,
  which = "all_genes",
  plot_type = "ranked_points",
  top_n = 50L,
  gene_column = NULL,
  score_column = NULL,
  highlight_genes = NULL,
  title = "conseguiR Diffusion Scores",
  subtitle = NULL,
  width = 10,
  height = 7,
  dpi = 300
) {
  if (missing(file_path) || is.null(file_path) || !nzchar(file_path)) {
    stop("`file_path` must be provided to save the diffusion plot.")
  }

  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)

  plot_obj <- create_diffusion_plot(
    diffusion = diffusion,
    table = table,
    which = which,
    plot_type = plot_type,
    top_n = top_n,
    gene_column = gene_column,
    score_column = score_column,
    highlight_genes = highlight_genes,
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
