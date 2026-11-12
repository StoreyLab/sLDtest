test_that("compute_entropy calculates valid Shannon entropy for standard inputs", {
  set.seed(42)
  p_vals <- runif(100)
  
  ent <- compute_entropy(p_vals, n_bin = 10, include_first_bin = TRUE)
  
  expect_true(is.numeric(ent))
  expect_length(ent, 1)
  expect_true(is.finite(ent))
  expect_true(ent >= 0) # Shannon entropy should be non-negative
})

test_that("compute_entropy respects the include_first_bin flag", {
  set.seed(123)
  # Create p-values where many are in the first bin ([0, 0.1])
  p_vals <- c(runif(80, 0, 0.1), runif(20, 0.1, 1))
  
  ent_with_first <- compute_entropy(p_vals, n_bin = 10, include_first_bin = TRUE)
  ent_without_first <- compute_entropy(p_vals, n_bin = 10, include_first_bin = FALSE)
  
  # Excluding the first bin should change the computed entropy value
  expect_true(ent_with_first != ent_without_first)
})

test_that("compute_entropy correctly validates input arguments", {
  p_vals <- runif(50)
  
  # Error when non-numeric p_vals
  expect_error(compute_entropy(as.character(p_vals), n_bin = 10),
               "'p_vals' must be a numeric vector."
               )
  
  # Error when p_vals outside the [0, 1] interval
  expect_error(compute_entropy(c(-0.01, p_vals[-1]), n_bin = 10),
               "All elements in 'p_vals' must lie between 0 and 1."
               )
  expect_error(compute_entropy(c(1.05, p_vals[-1]), n_bin = 10),
               "All elements in 'p_vals' must lie between 0 and 1."
               )
  
  # Error when invalid n_bin 
  expect_error(compute_entropy(p_vals, n_bin = 0),
               "'n_bin' must be a positive integer."
               )
  expect_error(compute_entropy(p_vals, n_bin = -3),
               "'n_bin' must be a positive integer."
               )
  expect_error(compute_entropy(p_vals, n_bin = 100), # 100 > length(p_vals) 
               "'n_bin' cannot exceed the number of elements in 'p_vals'."
               )
  expect_error(compute_entropy(p_vals, n_bin = "five"),
               "'n_bin' must be a positive integer."
               )
})