test_that("mu returns correctly", {
  result <- myncurve(mu = 10, sigma = 5, a = 6)
  expect_equal(result$mu, 10)
})

test_that("sigma returns correctly", {
  result <- myncurve(mu = 10, sigma = 5, a = 6)
  expect_equal(result$sigma, 5)
})

test_that("area returns correctly", {
  result <- myncurve(mu = 10, sigma = 5, a = 6)
  expect_equal(result$area, round(pnorm(6, 10, 5), 4))
})
