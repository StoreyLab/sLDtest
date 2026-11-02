test_that("sim_null_alt_pairs outputs correct dimensions and data types", {

  set.seed(42)
  I_loci <- 100
  J_ind <- 30
  Pi <- matrix(runif(I_loci * J_ind, min = 0.1, max = 0.9), nrow = I_loci, ncol = J_ind)
  
  # Choose indices for 5 alternate hypotheses
  alt_ids <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ncol = 2, byrow = TRUE)
  rho <- 0.5
  
  # Run function
  X <- sim_null_alt_pairs(Pi, alt_ids, rho)
  
  # Checks
  expect_true(is.matrix(X))
  expect_equal(nrow(X), I_loci)
  expect_equal(ncol(X), J_ind)
  expect_true(is.numeric(X))
  
  # Genotypes for biallelic markers should be counts (0, 1, 2)
  expect_true(all(X %in% c(0, 1, 2)))
})

test_that("sim_null_alt_pairs correctly validates input arguments", {
  Pi <- matrix(runif(50 * 10), nrow = 50, ncol = 10)
  alt_ids <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  rho <- 0.2
  
  # Error when Pi is not a matrix
  expect_error(
    sim_null_alt_pairs(Pi = data.frame(Pi), alt_ids = alt_ids, rho = rho),
    "'Pi' must be a numeric matrix."
  )
  
  # Error when alt_ids is not a matrix or doesn't have 2 columns
  expect_error(
    sim_null_alt_pairs(Pi = Pi, alt_ids = c(1, 2, 3), rho = rho),
    "'alt_ids' must be a two-column matrix of locus pairs."
  )
  expect_error(
    sim_null_alt_pairs(Pi = Pi, alt_ids = matrix(1:3, ncol = 3), rho = rho),
    "'alt_ids' must be a two-column matrix of locus pairs."
  )
  
  # Error when alt_ids indices exceed rows in Pi
  bad_alt_ids <- matrix(c(1, 2, 100, 200), ncol = 2, byrow = TRUE) # 100 & 200 > 50 rows
  expect_error(
    sim_null_alt_pairs(Pi = Pi, alt_ids = bad_alt_ids, rho = rho),
    "Indices in 'alt_ids' exceed the number of rows in 'Pi'."
  )
  
  # Error when rho is out of bounds or not numeric
  expect_error(
    sim_null_alt_pairs(Pi = Pi, alt_ids = alt_ids, rho = 1.5),
    "'rho' must be a single numeric value between -1 and 1."
  )
  expect_error(
    sim_null_alt_pairs(Pi = Pi, alt_ids = alt_ids, rho = "high"),
    "'rho' must be a single numeric value between -1 and 1."
  )
})