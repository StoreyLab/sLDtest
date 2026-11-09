#' sLD_SNP
#'
#' This function performs the sLD test for a given pair of SNPs at distinct loci \eqn{a} and \eqn{b},
#' returning the sLD test statistic. 
#'
#' @param X_a Length-\eqn{J} vector of genotypes at locus a. 
#' @param X_b Length-\eqn{J} vector of genotypes at locus b.
#' @param pi_a Length-\eqn{J} vector of IAFs at locus a.
#' @param pi_b Length-\eqn{J} vector of IAFs at locus b. 
#'
#' @return sLD test statistic computed for SNPs at loci a and b.
#'     
#' @note \code{X_a}, \code{X_b}, \code{pi_a}, and \code{pi_b} must all be of the same length.
#'
#' @examples
#' pi_a <- c(rep(0.25, 50), rep(0.3, 50))
#' pi_b <- c(rep(0.45, 50), rep(0.1, 50))
#' X_a <- rbinom(100, 2, pi_a)
#' X_b <- rbinom(100, 2, pi_b)
#' sLD_test_stat <- sLD_SNP(X_a, X_b, pi_a, pi_b) 
#'
#' @export

sLD_SNP <- function(X_a, X_b, pi_a, pi_b) {
  # Check arguments
  if (!is.numeric(X_a) || !is.numeric(X_b) || !is.numeric(pi_a) || !is.numeric(pi_b)) {
    stop("All input arguments must be numeric vectors.")
  }
  if (length(X_a) != length(X_b) || length(X_a) != length(pi_a) || length(X_a) != length(pi_b)) {
    stop("'X_a', 'X_b', 'pi_a', and 'pi_b' must have the same length.")
  }
  
  keep_id <- which(!is.na(X_a) & !is.na(X_b)) # check for missing genotype values
  n_ind <- length(keep_id)
  
  if (n_ind == 0) {
    warning("No valid non-missing pairs found. Returning NA.")
    return(NA_real_)
  }
  
  X_a <- X_a[keep_id]
  X_b <- X_b[keep_id]
  
  if (any(!X_a %in% c(0, 1, 2)) || any(!X_b %in% c(0, 1, 2))) {
    warning("Genotype values in 'X_a' or 'X_b' contain values outside {0, 1, 2}.")
  }
  
  pi_a <- pi_a[keep_id]
  pi_b <- pi_b[keep_id]
  
  mean_term <- (X_a - 2*pi_a) * (X_b - 2*pi_b)
  var_term <- 4*pi_a*(1 - pi_a)*pi_b*(1 - pi_b)
  
  denom <- sum(var_term)
  if (denom == 0) {
    warning("Variance term is zero. Returning NA.")
    return(NA_real_)
  }
  
  test_stat <- sum(mean_term) / sqrt(denom)
  
  return(test_stat)
}