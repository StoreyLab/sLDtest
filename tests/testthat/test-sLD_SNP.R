test_that("sLD_SNP calculates test statistic correctly for valid inputs", {
  set.seed(42)
  J <- 50
  pi_a <- runif(J, min = 0.1, max = 0.9)
  pi_b <- runif(J, min = 0.1, max = 0.9)
  X_a <- rbinom(J, size = 2, prob = pi_a)
  X_b <- rbinom(J, size = 2, prob = pi_b)
  
  stat <- sLD_SNP(X_a, X_b, pi_a, pi_b)
  
  expect_true(is.numeric(stat))
  expect_length(stat, 1)
  expect_true(is.finite(stat))
})

test_that("sLD_SNP handles missing (NA) genotype values correctly", {
  set.seed(123)
  J <- 30
  pi_a <- runif(J, min = 0.1, max = 0.9)
  pi_b <- runif(J, min = 0.1, max = 0.9)
  X_a <- rbinom(J, size = 2, prob = pi_a)
  X_b <- rbinom(J, size = 2, prob = pi_b)
  
  # Inject NAs
  X_a[c(3, 12)] <- NA
  X_b[c(12, 20)] <- NA
  
  stat <- sLD_SNP(X_a, X_b, pi_a, pi_b)
  
  expect_true(is.numeric(stat))
  expect_length(stat, 1)
  expect_true(is.finite(stat))
})

test_that("sLD_SNP handles zero denominator and zero valid individuals with warnings and NAs", {
  pi_a <- rep(0.5, 10)
  pi_b <- rep(0.5, 10)
  X_a <- rep(1, 10)
  X_b <- rep(1, 10)
  
  # Force zero variance denominator (e.g. all allele frequencies at 0 or 1)
  pi_a_zero <- rep(0, 10)
  expect_warning(stat_zero <- sLD_SNP(X_a, X_b, pi_a_zero, pi_b),
                 "Variance term is zero. Returning NA."
                 )
  expect_identical(stat_zero, NA_real_)
  
  # Force zero valid individuals (all NA)
  X_a_na <- rep(NA_real_, 10)
  expect_warning(stat_none <- sLD_SNP(X_a_na, X_b, pi_a, pi_b),
                 "No valid non-missing pairs found. Returning NA."
                 )
  expect_identical(stat_none, NA_real_)
})

test_that("sLD_SNP warns when genotypes are outside {0, 1, 2}", {
  pi_a <- rep(0.5, 10)
  pi_b <- rep(0.5, 10)
  X_a <- c(3, rep(1, 9)) # 3 is invalid for diploid genotype counts
  X_b <- rep(1, 10)
  
  expect_warning(sLD_SNP(X_a, X_b, pi_a, pi_b),
                 "Genotype values in 'X_a' or 'X_b' contain values outside \\{0, 1, 2\\}."
                 )
})

test_that("sLD_SNP correctly validates input arguments", {
  pi_a <- rep(0.5, 10)
  pi_b <- rep(0.5, 10)
  X_a <- rep(1, 10)
  X_b <- rep(1, 10)
  
  # Error when non-numeric input 
  expect_error(sLD_SNP(as.character(X_a), X_b, pi_a, pi_b),
               "All input arguments must be numeric vectors."
               )
  
  # Error when there is a length mismatch 
  expect_error(sLD_SNP(X_a, X_b[1:5], pi_a, pi_b),
               "'X_a', 'X_b', 'pi_a', and 'pi_b' must have the same length."
               )
})