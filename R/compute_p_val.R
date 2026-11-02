#' compute_p_val
#'
#' This function computes a p-value for a given sLD test statistic
#' assuming a theoretical Normal(0,1) null distribution. 
#' 
#' @param test_stats A vector of sLD test statistics.  
#' 
#' @return \eqn{p}-values corresponding to the given test statistic.  
#'     
#' @note Computes \eqn{p}-values using the theoretical Normal(0,1) CDF. 
#' 
#' @examples
#' pi_a <- c(rep(0.25, 50), rep(0.3, 50))
#' pi_b <- c(rep(0.45, 50), rep(0.1, 50))
#' X_a <- rbinom(100, 2, pi_a)
#' X_b <- rbinom(100, 2, pi_b)
#' sLD_test_stat <- sLD_SNP(X_a, X_b, pi_a, pi_b) 
#' sLD_p_val <- compute_p_val(sLD_test_stat)
#' 
#' @export

compute_p_val <- function(test_stats) {
  if (!is.numeric(test_stats)) {
    stop("'test_stats' must be a numeric vector.")
  }
  
  return(2 * pnorm(abs(test_stats), lower.tail = FALSE))
}