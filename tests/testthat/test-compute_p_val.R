test_that("compute_p_val calculates correct p-values for standard inputs", {
  # Test statistic of 0 should yield a p-value of 1
  expect_equal(compute_p_val(0), 1)
  
  # Test statistic of 1.96 should yield a p-value of about 0.05 (Normal CDF)
  expected_p_05 <- 2 * pnorm(-1.96)
  expect_equal(compute_p_val(1.96), expected_p_05, tolerance = 1e-4)
  
  # Negative test statistics should give the same p-values due to abs()
  expect_equal(compute_p_val(-1.96), compute_p_val(1.96))
  
  # Vector input should return a vector of identical length with values in [0, 1]
  stats <- c(0, 1.0, -1.96, 3.2, 5.0)
  p_vals <- compute_p_val(stats)
  
  expect_true(is.numeric(p_vals))
  expect_equal(length(p_vals), length(stats))
  expect_true(all(p_vals >= 0 & p_vals <= 1))
})

test_that("compute_p_val correctly checks input arguments", {
  # Error when non-numeric input 
  expect_error(compute_p_val("1.96"),
               "'test_stats' must be a numeric vector."
               )
  expect_error(compute_p_val(list(1, 2, 3)),
               "'test_stats' must be a numeric vector."
               )
})