#' compute_p_val_empr
#'
#' This function computes a two-sided p-value for each sLD test statistic
#' using an empirical null distribution
#' 
#' @param test_stats Vector of sLD test statistics. 
#' @param null_stats Vector of null test statistics. 
#' 
#' @return \eqn{p}-values corresponding to the given test statistics.  
#'     
#' @note Computes \eqn{p}-values using the empirical distribution corresponding to
#'      \code{null_stats}. 
#' 
#' @examples
#' pi_a <- c(rep(0.25, 50), rep(0.3, 50))
#' pi_b <- c(rep(0.45, 50), rep(0.1, 50))
#' X_a <- rbinom(100, 2, pi_a)
#' X_b <- rbinom(100, 2, pi_b)
#' test_stat <- sLD_SNP(X_a, X_b, pi_a, pi_b)
#' null_stats <- rnorm(1000) # arbitrary empirical null distribution
#' sLD_p_val <- compute_p_val_empr(test_stat, null_stats)
#' 
#' @export

compute_p_val_empr <- function(test_stats, null_stats) {
  # Check arguments 
  if (!is.numeric(test_stats) || !is.numeric(null_stats)) {
    stop("Both 'test_stats' and 'null_stats' must be numeric vectors.")
  }
  M <- length(null_stats)
  if (M == 0) {
    stop("The 'null_stats' vector cannot be empty.")
  }
  
  obs_eval <- abs(test_stats)
  null_eval <- abs(null_stats)
  null_sorted <- sort(null_eval)
  idx_below <- findInterval(obs_eval, null_sorted, left.open = TRUE)
  count_geq <- M - idx_below
  p_vals <- (count_geq + 1) / (M + 1) # add pseudocount 
  
  return(p_vals)
}