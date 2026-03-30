# Mapping scores to graph nodes

map_gene_scores_to_graph <- function(graph, gene_scores, ...) {
  stop("Not yet implemented: map_gene_scores_to_graph()")
}

map_reg_scores_to_graph <- function(graph, reg_scores, ...) {
  stop("Not yet implemented: map_reg_scores_to_graph()")
}

#' Assemble node prize table
#'
#' Combines modality-specific score tables into a unified node-level prize table.
#'
#' @param gene_scores Optional gene-level score table.
#' @param reg_scores Optional regulatory-element-level score table.
#' @param ... Additional layer tables.
#' @return A node prize table.
#' @export
assemble_node_prize_table <- function(gene_scores = NULL, reg_scores = NULL, ...) {
  stop("Not yet implemented: assemble_node_prize_table()")
}

fill_missing_node_scores <- function(x, fill = 0) {
  stop("Not yet implemented: fill_missing_node_scores()")
}
