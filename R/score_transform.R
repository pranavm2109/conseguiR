# Score transformation and normalization

transform_scores <- function(x, method = c("none", "clip", "sqrt", "asinh", "ecdf", "rank"), ...) {
  stop("Not yet implemented: transform_scores()")
}

normalize_scores_within_layer <- function(x, method = c("zscore", "rank", "minmax"), ...) {
  stop("Not yet implemented: normalize_scores_within_layer()")
}

standardize_scores_for_diffusion <- function(x, ...) {
  stop("Not yet implemented: standardize_scores_for_diffusion()")
}

combine_layer_scores <- function(...) {
  stop("Not yet implemented: combine_layer_scores()")
}
