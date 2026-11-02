#' sim_structural_LD
#'
#' This function draws pairs of genotypes \eqn{X_a}, \eqn{X_b} in structural LD
#' according to a multinomial distribution parameterized by individual-specific
#' haplotype frequencies (IHFs) with a targeted residual correlation \eqn{\rho}.
#' The residual correlation is applied between the minor alleles at loci \eqn{a}
#' and \eqn{b}.
#'
#' @param pi_a Length-\eqn{J} vector of individual-specific allele frequencies (IAF)
#'      at locus \eqn{a}.
#' @param pi_b Length-\eqn{J} vector of IAFs at locus \eqn{b}.
#' @param rho A single numeric value in \eqn{[-1, 1]} denoting the residual correlation
#'      between the minor alleles of loci \eqn{a} and \eqn{b}. Positive values indicate
#'      positive LD between minor alleles; negative values indicate negative LD.
#'
#' @return A \eqn{2 \times J} matrix where the first row contains genotypes for locus \eqn{a}
#'      and the second row contains genotypes for locus \eqn{b} across all \eqn{J} individuals,
#'      aligned with the input reference orientation of \code{pi_a} and \code{pi_b}.
#'
#' @note \code{pi_a} and \code{pi_b} must have the same length \eqn{J}. \code{rho} must
#'      be in \eqn{[-1, 1]}. 
#'
#' @examples
#' pi_a <- c(rep(0.75, 50), rep(0.8, 50)) 
#' pi_b <- c(rep(0.2, 50), rep(0.15, 50)) 
#' rho <- 0.2
#' X_ab <- sim_structural_LD(pi_a, pi_b, rho = rho)
#'
#' @export
sim_structural_LD <- function(pi_a, pi_b, rho) {
  # Check arguments
  if (!is.numeric(pi_a) || !is.numeric(pi_b)) {
    stop("'pi_a' and 'pi_b' must be numeric vectors.")
  }
  if (length(pi_a) != length(pi_b)) {
    stop("'pi_a' and 'pi_b' must be of the same length.")
  }
  if (any(pi_a < 0 | pi_a > 1) || any(pi_b < 0 | pi_b > 1)) {
    stop("Allele frequencies must lie between 0 and 1.")
  }
  if (!is.numeric(rho) || length(rho) != 1 || rho < -1 || rho > 1) {
    stop("'rho' must be a single numeric value between -1 and 1.")
  }
  
  # Check whether each locus tracks the major allele; if so, flip to minor allele
  flip_a <- mean(pi_a) > 0.5
  flip_b <- mean(pi_b) > 0.5
  
  p_a_eff <- if (flip_a) 1 - pi_a else pi_a
  p_b_eff <- if (flip_b) 1 - pi_b else pi_b
  
  # Compute individual-specific LD param
  sd_prod <- sqrt(p_a_eff * (1 - p_a_eff) * p_b_eff * (1 - p_b_eff))
  D_target <- rho * sd_prod
  
  # Compute theoretical bounds on LD param
  D_lower <- pmax(-p_a_eff * p_b_eff, -(1 - p_a_eff) * (1 - p_b_eff))
  D_upper <- pmin(p_a_eff * (1 - p_b_eff), (1 - p_a_eff) * p_b_eff)
  
  # Truncate LD param if necessary
  eps <- 1e-12
  LD_offset <- pmin(pmax(D_target, D_lower + eps), D_upper - eps)
  
  # Compute haplotypes
  h_11 <- pmax(0, p_a_eff * p_b_eff + LD_offset)
  h_10 <- pmax(0, p_a_eff * (1 - p_b_eff) - LD_offset)
  h_01 <- pmax(0, (1 - p_a_eff) * p_b_eff - LD_offset)
  h_00 <- pmax(0, (1 - p_a_eff) * (1 - p_b_eff) + LD_offset)
  prob_sums <- h_11 + h_10 + h_01 + h_00
  haplo_freqs <- cbind(h_11 / prob_sums,
                       h_10 / prob_sums,
                       h_01 / prob_sums,
                       h_00 / prob_sums)
  
  # Draw haplotypes and genotypes
  n_ind <- length(pi_a)
  haplo <- t(sapply(1:n_ind, function(i) rmultinom(1, size = 2, prob = haplo_freqs[i, ])))
  X_a_eff <- haplo[, 1] + haplo[, 2]
  X_b_eff <- haplo[, 1] + haplo[, 3]
  X_ab <- rbind(X_a_eff, X_b_eff)
  
  return(X_ab)
}