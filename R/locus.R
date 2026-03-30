# Locus-level interpretation

rank_genes_by_diffusion_gain <- function(x, ...) {
  stop("Not yet implemented: rank_genes_by_diffusion_gain()")
}

#' Identify standout regulatory elements
#'
#' Identifies regulatory elements that stand out for a gene relative to its other linked elements.
#'
#' @param gene_reg_scores A gene-regulatory score table.
#' @param ... Additional arguments controlling standout selection.
#' @return A table of standout regulatory elements.
#' @export
identify_standout_reg_elements <- function(gene_reg_scores, ...) {
  stop("Not yet implemented: identify_standout_reg_elements()")
}

summarize_gene_locus_support <- function(gene, scores, ...) {
  stop("Not yet implemented: summarize_gene_locus_support()")
}

#' Extract a gene-regulatory neighborhood
#'
#' Extracts the local graph neighborhood around a selected gene and its linked regulatory elements.
#'
#' @param graph A graph object.
#' @param gene A gene identifier.
#' @param ... Additional arguments controlling the neighborhood extraction.
#' @return A neighborhood subgraph or annotated table.
#' @export
extract_gene_reg_neighborhood <- function(graph, gene, ...) {
  stop("Not yet implemented: extract_gene_reg_neighborhood()")
}

#' Run a locus report pipeline
#'
#' Generates a locus-level report with linked elements, standout enhancers, and track-ready outputs.
#'
#' @param gene A gene identifier.
#' @param ... Additional inputs required by the implementation.
#' @return A structured locus report object.
#' @export
run_locus_report_pipeline <- function(gene, ...) {
  stop("Not yet implemented: run_locus_report_pipeline()")
}
