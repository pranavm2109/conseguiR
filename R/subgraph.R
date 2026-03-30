# Subgraph calling

#' Call a cardinality-constrained subgraph
#'
#' Selects a fixed-size subgraph that balances node prize and edge quality.
#'
#' @param graph A graph object.
#' @param node_prizes Node prize table.
#' @param k Number of nodes to select.
#' @param ... Additional arguments controlling the optimizer.
#' @return A selected subgraph object or list.
#' @export
call_subgraph_cardinality <- function(graph, node_prizes, k = 50, ...) {
  stop("Not yet implemented: call_subgraph_cardinality()")
}

#' Call a prize-collecting Steiner subgraph
#'
#' Selects a connected subgraph or forest using a PCST-style objective.
#'
#' @param graph A graph object.
#' @param node_prizes Node prize table.
#' @param ... Additional arguments controlling the solver.
#' @return A selected subgraph object or list.
#' @export
call_subgraph_pcst <- function(graph, node_prizes, ...) {
  stop("Not yet implemented: call_subgraph_pcst()")
}

call_subgraph_forest <- function(graph, node_prizes, ...) {
  stop("Not yet implemented: call_subgraph_forest()")
}

rank_candidate_subgraphs <- function(subgraphs, ...) {
  stop("Not yet implemented: rank_candidate_subgraphs()")
}

compare_subgraph_solutions <- function(...) {
  stop("Not yet implemented: compare_subgraph_solutions()")
}
