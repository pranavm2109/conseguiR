# Somatic scoring

#' Score somatic signal at genes with dndscv
#'
#' Computes gene-level somatic scores from a somatic mutation table using dndscv.
#'
#' @param maf Somatic mutation table.
#' @param ... Additional arguments forwarded to the implementation.
#' @return A gene-level somatic score table.
#' @export
score_somatic_genes_dndscv <- function(maf, ...) {
  stop("Not yet implemented: score_somatic_genes_dndscv()")
}

#' Score somatic signal at regulatory elements with fishHook
#'
#' Computes regulatory-element-level somatic enrichment scores using fishHook.
#'
#' @param maf Somatic mutation table.
#' @param reg_elements Regulatory element table.
#' @param covariates Optional covariate table.
#' @param ... Additional arguments forwarded to the implementation.
#' @return A regulatory-element-level somatic score table.
#' @export
score_somatic_regs_fishhook <- function(maf, reg_elements, covariates = NULL, ...) {
  stop("Not yet implemented: score_somatic_regs_fishhook()")
}

prepare_fishhook_covariates <- function(reg_elements, aid = NULL, ...) {
  stop("Not yet implemented: prepare_fishhook_covariates()")
}

score_somatic_regs_fishhook_aid <- function(maf, reg_elements, aid, ...) {
  stop("Not yet implemented: score_somatic_regs_fishhook_aid()")
}
