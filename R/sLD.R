#' sLD
#'
#' This function performs the sLD test for a set of SNP pairs, returning a 
#' corresponding set of test statistics. 
#'
#' @param X Genotype matrix (\eqn{I} loci x \eqn{J} individuals).
#' @param Pi \eqn{I} x \eqn{J} matrix of individual-specific allele frequencies. 
#' @param pair_ids Pairs of locus indices to test. 
#'
#' @return sLD test statistics for all locus pairs in pair_ids. 
#'     
#' @note X and Pi must have the same dimensions. 
#'      pair_ids should be a matrix of dimension # pairs x 2. 
#'
#' @examples
#' Pi <- matrix(runif(n = 500 * 50), nrow = 500, ncol = 50)
#' X <- matrix(rbinom(500 * 50, size = 2, prob = Pi), nrow = 500, ncol = 50)
#' pair_ids <- t(replicate(10, sample(1:500, size = 2, replace = FALSE)))
#' sLD_test_stat <- sLD(X, Pi, pair_ids) 
#' 
#' @export

sLD <- function(X, Pi, pair_ids) {
  # Check arguments
  if (!is.matrix(X) || !is.matrix(Pi)) {
    stop("'X' and 'Pi' must be numeric matrices.")
  }
  if (!any(dim(X) == dim(Pi))) {
    stop("'X' and 'Pi' must be of the same size.")
  }
  if (!is.matrix(pair_ids) || ncol(pair_ids) != 2) {
    stop("'pair_ids' must be a two-column matrix of locus pairs.")
  }
  if (max(pair_ids) > nrow(X) || min(pair_ids) < 1) {
    stop("Indices in 'pair_ids' exceed the row dimensions of 'X' and 'Pi'.")
  }
  
  n_pairs <- nrow(pair_ids)
  test_stat_vec <- numeric(n_pairs)
  
  for (pair in 1:n_pairs) {
    locus_a <- pair_ids[pair, 1]
    locus_b <- pair_ids[pair, 2]
    X_a <- X[locus_a, ]
    X_b <- X[locus_b, ]
    pi_a <- Pi[locus_a, ]
    pi_b <- Pi[locus_b, ]
    test_stat_vec[pair] <- sLD_SNP(X_a, X_b, pi_a, pi_b) 
  }
  
  return(test_stat_vec)
}