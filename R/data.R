#' Simulated Genotype Matrix
#'
#' A simulated genotype matrix containing individuals and loci used as 
#' an example dataset for the `sLDtest` package. The genotypes were simulated 
#' under an admixture model with $K=4$ admixed populations.
#'
#' @format Numeric matrix with $I$ rows (loci) and $J$ columns (individuals). 
#'   Values represent counts of the allele coded as "1" at each locus ($0, 1, 2$).
#' @source Simulated using the `bnpsd` package.
#' @usage data(genodata)
#' @keywords datasets
"genodata"