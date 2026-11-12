test_that("compute_p_val_empr calculates correct p-values when there are ties and duplicates", {
  set.seed(42)
  
  # Sample null stats 
  null_stats <- c(-4, -2, 0, 1, 3) # M = 5 null stats
  
  # What if there is an exact match to an absolute null value
  p_val_exact <- compute_p_val_empr(2, null_stats)
  expect_equal(p_val_exact, 4 / 6)
  
  # Check that the function correctly handles negative test statistics
  p_val_neg <- compute_p_val_empr(-2, null_stats)
  expect_equal(p_val_neg, p_val_exact)
  
  # Check inputs and outputs when inputting vector of observed test stats
  test_stats <- c(0, 5, -3)
  p_vals <- compute_p_val_empr(test_stats, null_stats)
  expect_true(is.numeric(p_vals))
  expect_equal(length(p_vals), length(test_stats))
  expect_equal(p_vals, c(1.0, 1 / 6, 0.5))
  expect_true(all(p_vals >= 0 & p_vals <= 1))
  
  # Handling duplicates in null_stats:
  dup_nulls <- c(-2, 2, -2, 0, 0)
  
  # Test obs = 2: should be 3 null values >= 2
  # p-value = (3 + 1) / (5 + 1) = 4 / 6 = 2/3
  expect_equal(compute_p_val_empr(2, dup_nulls), 2/3)
  
  # Test obs = 0: should be 5 null values >= 0 
  # p-value = (5 + 1) / (5 + 1) = 1
  expect_equal(compute_p_val_empr(0, dup_nulls), 1)
  
  # Passing identical observed test statistics should result in identical p-values
  repeated_tests <- c(2, 2, 2, -2)
  repeated_p_vals <- compute_p_val_empr(repeated_tests, null_stats)
  expect_equal(length(repeated_p_vals), 4)
  expect_equal(unique(repeated_p_vals), 4 / 6)
  
  # Extreme values outside the null distribution range
  expect_equal(compute_p_val_empr(100, null_stats), 1 / (length(null_stats) + 1))
})

test_that("compute_p_val_empr correctly validates input arguments", {
  null_stats <- c(1, 2, 3)
  
  # Error when non-numeric test_stats or null_stats
  expect_error(compute_p_val_empr("1.5", null_stats),
               "Both 'test_stats' and 'null_stats' must be numeric vectors."
               )
  expect_error(compute_p_val_empr(1.5, "not_numeric"),
               "Both 'test_stats' and 'null_stats' must be numeric vectors."
              )
  
  # Error when empty null_stats vector
  expect_error(compute_p_val_empr(1.5, numeric(0)),
               "The 'null_stats' vector cannot be empty."
               )
  })