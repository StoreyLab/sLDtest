test_that("sim_structural_LD outputs correct dimensions and valid genotype values", {
  set.seed(123)
  J_ind <- 50
  pi_a <- runif(J_ind, min = 0.1, max = 0.9)
  pi_b <- runif(J_ind, min = 0.1, max = 0.9)
  rho <- 0.3
  
  X_ab <- sim_structural_LD(pi_a, pi_b, rho)
  
  # Checks
  expect_true(is.matrix(X_ab))
  expect_equal(nrow(X_ab), 2)
  expect_equal(ncol(X_ab), J_ind)
  expect_true(is.numeric(X_ab))
  expect_true(all(is.finite(X_ab)))
  
  # Genotypes for biallelic markers should be counts (0, 1, 2)
  expect_true(all(X_ab %in% c(0, 1, 2)))
})

test_that("sim_structural_LD correctly validates input arguments", {
  pi_a <- rep(0.5, 10)
  pi_b <- rep(0.4, 10)
  rho <- 0.2
  
  # Error when non-numeric vector inputs
  expect_error(
    sim_structural_LD(as.character(pi_a), pi_b, rho),
    "'pi_a' and 'pi_b' must be numeric vectors."
  )
  expect_error(
    sim_structural_LD(pi_a, as.character(pi_b), rho),
    "'pi_a' and 'pi_b' must be numeric vectors."
  )
  
  # Error when there is a ength mismatch
  expect_error(
    sim_structural_LD(pi_a, pi_b[1:5], rho),
    "'pi_a' and 'pi_b' must be of the same length."
  )
  
  # Error when allele frequencies out of bounds [0, 1]
  expect_error(
    sim_structural_LD(c(-0.01, pi_a[-1]), pi_b, rho),
    "Allele frequencies must lie between 0 and 1."
  )
  expect_error(
    sim_structural_LD(pi_a, c(pi_b[-10], 1.05), rho),
    "Allele frequencies must lie between 0 and 1."
  )
  
  # Error when invalid rho
  expect_error(
    sim_structural_LD(pi_a, pi_b, 1.2),
    "'rho' must be a single numeric value between -1 and 1."
  )
  expect_error(
    sim_structural_LD(pi_a, pi_b, -1.5),
    "'rho' must be a single numeric value between -1 and 1."
  )
  expect_error(
    sim_structural_LD(pi_a, pi_b, "invalid_rho"),
    "'rho' must be a single numeric value between -1 and 1."
  )
  expect_error(
    sim_structural_LD(pi_a, pi_b, c(0.1, 0.5)),
    "'rho' must be a single numeric value between -1 and 1."
  )
})