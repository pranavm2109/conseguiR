# Epigenomic scoring

#' Score epigenomic signal at regulatory elements
#'
#' Computes regulatory-element-level epigenomic scores from one or more tracks.
#'
#' @param tracks Epigenomic tracks.
#' @param reg_elements Regulatory element table.
#' @param ... Additional arguments forwarded to the implementation.
#' @return A regulatory-element-level epigenomic score table.
#' @export
score_epigenomic_regs <- function(tracks, reg_elements, ...) {
  stop("Not yet implemented: score_epigenomic_regs()")
}

score_epigenomic_overlap <- function(tracks, reg_elements, ...) {
  stop("Not yet implemented: score_epigenomic_overlap()")
}

score_epigenomic_variability <- function(tracks, reg_elements, ...) {
  stop("Not yet implemented: score_epigenomic_variability()")
}

aggregate_epigenomic_tracks <- function(tracks, ...) {
  stop("Not yet implemented: aggregate_epigenomic_tracks()")
}
