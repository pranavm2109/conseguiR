# Plotting and visualization

plot_gene_scores <- function(x, ...) {
  stop("Not yet implemented: plot_gene_scores()")
}

plot_reg_scores <- function(x, ...) {
  stop("Not yet implemented: plot_reg_scores()")
}

#' Plot diffusion shift
#'
#' Visualizes how scores or ranks changed before and after diffusion.
#'
#' @param x Diffusion results table.
#' @param ... Additional plotting arguments.
#' @return A ggplot object.
#' @export
plot_diffusion_shift <- function(x, ...) {
  stop("Not yet implemented: plot_diffusion_shift()")
}

#' Plot a selected subgraph
#'
#' Draws a selected graph or subgraph using node and edge attributes.
#'
#' @param graph A graph object.
#' @param ... Additional plotting arguments.
#' @return A plot object.
#' @export
plot_subgraph <- function(graph, ...) {
  stop("Not yet implemented: plot_subgraph()")
}

#' Plot locus tracks
#'
#' Plots locus-level multimodal evidence across a genomic region.
#'
#' @param locus_data A locus-level data object.
#' @param ... Additional plotting arguments.
#' @return A plot object.
#' @export
plot_locus_tracks <- function(locus_data, ...) {
  stop("Not yet implemented: plot_locus_tracks()")
}

plot_modality_convergence <- function(x, ...) {
  stop("Not yet implemented: plot_modality_convergence()")
}

#' Plot expression by regulatory mutation status
#'
#' Compares expression between samples with mutations in selected regulatory elements and wild type samples.
#'
#' @param expression_data Expression matrix or table.
#' @param mutation_labels Mutation status labels.
#' @param ... Additional plotting arguments.
#' @return A ggplot object.
#' @export
plot_expression_by_re_mutation <- function(expression_data, mutation_labels, ...) {
  stop("Not yet implemented: plot_expression_by_re_mutation()")
}
