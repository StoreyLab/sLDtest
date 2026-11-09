#' compute_entropy
#'
#' This function computes the Shannon entropy of a set of \eqn{p}-values,
#' approximating how close their distribution is to the Uniform(0,1)
#' distribution. 
#' 
#' @param p_vals A vector of \eqn{p}-values.
#' @param n_bin Integer, the number of bins for the histogram.
#' @param include_first_bin Boolean, denotes whether to compute entropy 
#'   using the first interval (closest to 0) or not. 
#' 
#' @return Shannon entropy of the \eqn{p}-values.
#' 
#' @note \code{include_first_bin} is TRUE by default. 
#' 
#' @examples
#' p_vals <- runif(1000)
#' compute_entropy(p_vals, n_bin = 10, include_first_bin = TRUE)
#' 
#' @importFrom graphics hist
#' @export

compute_entropy <- function(p_vals, n_bin, include_first_bin = TRUE) {
  # Check arguments
  if (!is.numeric(p_vals)) {
    stop("'p_vals' must be a numeric vector.")
  }
  if (any(p_vals < 0 | p_vals > 1, na.rm = TRUE)) {
    stop("All elements in 'p_vals' must lie between 0 and 1.")
  }
  if (!is.numeric(n_bin) || length(n_bin) != 1 || n_bin <= 0) {
    stop("'n_bin' must be a positive integer.")
  }
  if (n_bin > length(p_vals)) {
    stop("'n_bin' cannot exceed the number of elements in 'p_vals'.")
  }
  
  h <- hist(p_vals, breaks = seq(0, 1, length.out = n_bin + 1), plot = FALSE)
  f <- h$counts / sum(h$counts) # get proportion in each bin
  
  entropy <- 0
  start_idx <- if (include_first_bin) 1 else 2
  
  for (i in start_idx:n_bin) {
    if (f[i] > 0) {
      entropy <- entropy + f[i] * log2(f[i])
    }
  }
  
  return(-entropy)
}