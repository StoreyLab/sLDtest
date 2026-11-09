#' Simulated Genotype Matrix
#'
#' A simulated genotype matrix containing individuals and loci used as 
#' an example dataset for the `sLDtest` package. The genotypes were simulated 
#' under an admixture model with \eqn{K=4} admixed populations.
#'
#' @format Numeric matrix with \eqn{I} rows (loci) and \eqn{J} columns (individuals). 
#'   Values represent counts of the allele coded as "1" at each locus (\eqn{0, 1, 2}).
#' @source Simulated using the `bnpsd` package.
#' @usage data(genodata)
#' @keywords datasets
"genodata"