# Graph construction and annotation

#' Build a gene-gene graph
#'
#' Constructs a gene-gene graph from STRING or another interaction backend.
#'
#' @param interactions Interaction table.
#' @return A graph object.
#' @export
build_gene_gene_graph <- function(interactions) {
  stop("Not yet implemented: build_gene_gene_graph()")
}

#' Build a gene-regulatory graph
#'
#' Constructs a bipartite graph linking genes and regulatory elements.
#'
#' @param links Gene-regulatory interaction table.
#' @return A graph object.
#' @export
build_gene_reg_graph <- function(links) {
  stop("Not yet implemented: build_gene_reg_graph()")
}

annotate_graph_nodes <- function(graph, node_metadata) {
  stop("Not yet implemented: annotate_graph_nodes()")
}

annotate_graph_edges <- function(graph, edge_metadata) {
  stop("Not yet implemented: annotate_graph_edges()")
}

subset_graph_by_nodes <- function(graph, nodes) {
  stop("Not yet implemented: subset_graph_by_nodes()")
}
