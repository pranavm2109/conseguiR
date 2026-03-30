# Germline scoring with MAGMA

#' Score germline signal at genes with MAGMA
#'
#' Computes gene-level germline scores from GWAS summary statistics using MAGMA.
#'
#' @param gwas GWAS summary statistics table.
#' @param annotation Optional precomputed MAGMA annotation.
#' @param ... Additional arguments forwarded to the implementation.
#' @return A gene-level germline score table.
#' @export
score_germline_genes_magma <- function(gwas, annotation = NULL, ...) {
  stop("Not yet implemented: score_germline_genes_magma()")
}

#' Score germline signal at regulatory elements with MAGMA
#'
#' Computes regulatory-element-level germline scores from GWAS summary statistics using MAGMA.
#'
#' @param gwas GWAS summary statistics table.
#' @param annotation Optional precomputed MAGMA annotation.
#' @param ... Additional arguments forwarded to the implementation.
#' @return A regulatory-element-level germline score table.
#' @export
score_germline_regs_magma <- function(gwas, annotation = NULL, ...) {
  stop("Not yet implemented: score_germline_regs_magma()")
}

prepare_magma_annotation_genes <- function(genes, snps, ...) {
  stop("Not yet implemented: prepare_magma_annotation_genes()")
}

prepare_magma_annotation_regs <- function(reg_elements, snps, ...) {
  stop("Not yet implemented: prepare_magma_annotation_regs()")
}

run_magma_gene_step <- function(...) {
  stop("Not yet implemented: run_magma_gene_step()")
}

run_magma_reg_step <- function(...) {
  stop("Not yet implemented: run_magma_reg_step()")
}
