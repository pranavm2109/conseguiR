# High-level wrappers

#' Run the full scoring pipeline
#'
#' Runs somatic, germline, and epigenomic scoring from raw inputs to standardized layer outputs.
#'
#' @param somatic Somatic input.
#' @param gwas Germline GWAS input.
#' @param epigenomic Epigenomic input.
#' @param ... Additional arguments controlling the pipeline.
#' @return A structured list of scored outputs.
#' @export
run_full_scoring_pipeline <- function(somatic, gwas, epigenomic, ...) {
  stop("Not yet implemented: run_full_scoring_pipeline()")
}

#' Run the graph integration pipeline
#'
#' Maps scores to graph nodes, runs diffusion, and returns integrated graph-level results.
#'
#' @param graph A graph object.
#' @param scores A scored input object.
#' @param ... Additional arguments controlling the pipeline.
#' @return A structured graph integration result.
#' @export
run_graph_integration_pipeline <- function(graph, scores, ...) {
  stop("Not yet implemented: run_graph_integration_pipeline()")
}

#' Run the subgraph pipeline
#'
#' Runs subgraph selection on integrated node prizes and graph edges.
#'
#' @param graph A graph object.
#' @param integrated_scores Integrated graph scores.
#' @param ... Additional arguments controlling the pipeline.
#' @return A selected subgraph object.
#' @export
run_subgraph_pipeline <- function(graph, integrated_scores, ...) {
  stop("Not yet implemented: run_subgraph_pipeline()")
}
