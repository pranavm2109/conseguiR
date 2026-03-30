# Diffusion and propagation

run_diffusion <- function(graph, node_scores, ...) {
  stop("Not yet implemented: run_diffusion()")
}

#' Run multilayer diffusion
#'
#' Performs graph diffusion separately or jointly across multiple score layers.
#'
#' @param graph A graph object.
#' @param node_prizes A unified node prize table.
#' @param ... Additional arguments controlling diffusion.
#' @return A table or object containing propagated scores.
#' @export
run_multilayer_diffusion <- function(graph, node_prizes, ...) {
  stop("Not yet implemented: run_multilayer_diffusion()")
}

transfer_reg_signal_to_genes <- function(graph, reg_scores, ...) {
  stop("Not yet implemented: transfer_reg_signal_to_genes()")
}

transfer_gene_signal_to_regs <- function(graph, gene_scores, ...) {
  stop("Not yet implemented: transfer_gene_signal_to_regs()")
}

inspect_diffusion_results <- function(x, ...) {
  stop("Not yet implemented: inspect_diffusion_results()")
}
