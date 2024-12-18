test_that("log_sum_exp works", {
  res1 <- log(exp(1) + exp(2) + exp(3))
  expect_equal(log_sum_exp(1, 2, 3), res1)
  expect_equal(log_sum_exp(c(1, 2, 3)), res1)
  # Large numbers should be finite and not compute exp directly
  expect_true(is.finite(log_sum_exp(1e6, 1e6 + 2, 1e6 + 10)))
})
