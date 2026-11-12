#' sim_null_alt_pairs
#'
#' This function draws a genotype matrix consisting of a mixture of genotype 
#' pairs in structural LD (true alternate hypotheses) and not in structural 
#' LD (true null hypotheses). 
#'
#' @param Pi Matrix of individual-specific allele frequencies 
#'      (\eqn{I} loci x \eqn{J} individuals).
#' @param alt_ids Pairs of locus indices corresponding to true alternate hypotheses. 
#' @param rho A single numeric value in \eqn{[-1, 1]} denoting the residual correlation
#'      between the minor alleles of loci \eqn{a} and \eqn{b}. Positive values indicate
#'      positive LD between minor alleles; negative values indicate negative LD.
#'
#' @return An \eqn{I} x \eqn{J} genotype matrix.
#'     
#' @note \code{alt_ids} should be a matrix of size # alternate hypotheses x \eqn{2}. 
#'      \code{rho} must be in the interval \eqn{[-1, 1]}. 
#'
#' @examples
#' Pi <- matrix(runif(n = 500*50), nrow = 500, ncol = 50)
#' alt_ids <- t(replicate(1000, sample(1:500, size = 2, replace = FALSE)))
#' rho <- 0.4 
#' X_ab <- sim_null_alt_pairs(Pi, alt_ids, rho)
#' 
#' @importFrom stats rbinom
#' @export

sim_null_alt_pairs <- function(Pi, alt_ids, rho) {
  # Check arguments
  if (!is.matrix(Pi)) {
    stop("'Pi' must be a numeric matrix.")
  }
  if (!is.matrix(alt_ids) || ncol(alt_ids) != 2) {
    stop("'alt_ids' must be a two-column matrix of locus pairs.")
  }
  
  I_loci <- nrow(Pi)
  J_ind <- ncol(Pi)
  
  if (max(alt_ids) > I_loci || min(alt_ids) < 1) {
    stop("Indices in 'alt_ids' exceed the number of rows in 'Pi'.")
  }
  if (!is.numeric(rho) || length(rho) != 1 || rho < -1 || rho > 1) {
    stop("'rho' must be a single numeric value between -1 and 1.")
  }
  
  n_pairs <- nrow(alt_ids)
  X <- matrix(nrow = I_loci, ncol = J_ind)
  non_alt_ids <- setdiff(1:I_loci, unique(as.matrix(alt_ids))) # get null indices
  
  # Draw genotype pairs in structural LD
  for (i in 1:n_pairs) {
    locus_a <- alt_ids[i, 1]
    locus_b <- alt_ids[i, 2]
    pi_a <- Pi[locus_a, ]
    pi_b <- Pi[locus_b, ]
    X[c(locus_a, locus_b), ] <- sim_structural_LD(pi_a, pi_b, rho)
  }
  
  # Draw genotypes at remaining loci (not in structural LD)
  for (j in non_alt_ids) {
    X[j, ] <- rbinom(J_ind, size = 2, prob = Pi[j, ])
  }
  
  return(X)
}