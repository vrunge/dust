### variance model: y close to 0 after large values
### (with plain cumulative sums, y^2 is lost and we get a false single point segment)
test_that("variance model: no false change for a value close to 0", {
  set.seed(1)
  y <- c(rnorm(10, sd = 1e6), rnorm(990))
  y[500] <- 1e-9
  penalty <- 4 * log(1000)   # the one-point segment would be optimal with 2 log(n)
  for (method in c("DUST", "DUSTib", "OP"))
    expect_equal(dust.1D(y, penalty, "variance", method)$changepoints, c(10, 1000),
                 info = method)
 
    expect_equal(dust.MD(rbind(y, y), 2 * penalty, "variance")$changepoints,
                 c(10, 1000))
  expect_true(is.finite(segmentation_Cost_1D(y, c(10, 499, 500, 1000), "variance")))
})
