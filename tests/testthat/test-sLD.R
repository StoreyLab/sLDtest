test_that("sLD calculates test statistics correctly for multiple pairs", {
  set.seed(42)
  I_loci <- 50
  J_ind <- 30
  
  Pi <- matrix(runif(I_loci * J_ind, min = 0.1, max = 0.9), nrow = I_loci, ncol = J_ind)
  X <- matrix(rbinom(I_loci * J_ind, size = 2, prob = Pi), nrow = I_loci, ncol = J_ind)
  
  # Select 5 pairs to test
  pair_ids <- matrix(c(1, 2,  5, 10,  12, 15,  20, 25,  30, 35), ncol = 2, byrow = TRUE)
  
  stats <- sLD(X, Pi, pair_ids)
  
  expect_true(is.numeric(stats))
  expect_length(stats, nrow(pair_ids))
  expect_true(all(is.finite(stats) | is.na(stats)))
})

test_that("sLD correctly validates input arguments", {
  I_loci <- 30
  J_ind <- 20
  Pi <- matrix(runif(I_loci * J_ind), nrow = I_loci, ncol = J_ind)
  X <- matrix(rbinom(I_loci * J_ind, size = 2, prob = Pi), nrow = I_loci, ncol = J_ind)
  pair_ids <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  
  # Error when non-matrix X or Pi 
  expect_error(
    sLD(data.frame(X), Pi, pair_ids),
    "'X' and 'Pi' must be numeric matrices."
  )
  expect_error(
    sLD(X, data.frame(Pi), pair_ids),
    "'X' and 'Pi' must be numeric matrices."
  )
  
  # Error when there is a dimension mismatch between X and Pi
  Pi_bad <- matrix(runif(40), nrow = 20, ncol = 2)
  expect_error(
    sLD(X, Pi_bad, pair_ids),
    "'X' and 'Pi' must be of the same size."
  )
  
  # Error when pair_ids is not of the correct format/type
  expect_error(
    sLD(X, Pi, c(1, 2, 3)),
    "'pair_ids' must be a two-column matrix of locus pairs."
  )
  expect_error(
    sLD(X, Pi, matrix(1:3, ncol = 3)),
    "'pair_ids' must be a two-column matrix of locus pairs."
  )
  
  # Error when indices exceed row dimensions or are out of bounds
  bad_pair_high <- matrix(c(1, 2, 5, I_loci + 1), ncol = 2, byrow = TRUE)
  expect_error(
    sLD(X, Pi, bad_pair_high),
    "Indices in 'pair_ids' exceed the row dimensions of 'X' and 'Pi'."
  )
  
  bad_pair_low <- matrix(c(0, 1, 2, 3), ncol = 2, byrow = TRUE)
  expect_error(
    sLD(X, Pi, bad_pair_low),
    "Indices in 'pair_ids' exceed the row dimensions of 'X' and 'Pi'."
  )
})