### variance model: y close to 0 after large values
### (with plain cumulative sums, y^2 is lost and we get a false single point segment)
test_that("variance model: no false change for a value close to 0", {
  set.seed(1)
  y <- c(rnorm(10, sd = 1e6), rnorm(990))
  y[500] <- 1e-9
  penalty <- 2 * log(1000)
  for (backend in c("scalar", "highway")) for (method in c("DUST", "DUSTib", "OP"))
    expect_equal(dust.1D(y, penalty, "variance", method, backend)$changepoints, c(10, 1000),
                 info = paste(method, backend))
  for (backend in c("scalar", "highway"))
    expect_equal(dust.MD(rbind(y, y), 2 * penalty, "variance", backend = backend)$changepoints,
                 c(10, 1000), info = backend)
  expect_true(is.finite(segmentation_Cost_1D(y, c(10, 499, 500, 1000), "variance")))
})
