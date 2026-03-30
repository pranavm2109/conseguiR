# Export and app integration

#' Export scores for app consumption
#'
#' Writes processed score tables in an app-friendly format.
#'
#' @param x A score object or table.
#' @param path Output path.
#' @return Invisibly returns the output path.
#' @export
export_scores_for_app <- function(x, path, ...) {
  stop("Not yet implemented: export_scores_for_app()")
}

#' Export a selected subgraph for app consumption
#'
#' Writes selected nodes and edges in an app-friendly format.
#'
#' @param x A subgraph object.
#' @param path Output path.
#' @return Invisibly returns the output path.
#' @export
export_subgraph_for_app <- function(x, path, ...) {
  stop("Not yet implemented: export_subgraph_for_app()")
}

#' Export locus tracks
#'
#' Writes track files such as BED and bedGraph for genome browser or app visualization.
#'
#' @param x A locus-level data object.
#' @param path Output directory.
#' @return Invisibly returns the output path.
#' @export
export_locus_tracks <- function(x, path, ...) {
  stop("Not yet implemented: export_locus_tracks()")
}

export_standout_enhancer_lists <- function(x, path, ...) {
  stop("Not yet implemented: export_standout_enhancer_lists()")
}
