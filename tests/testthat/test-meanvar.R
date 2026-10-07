test_that("meanVar methods and backends give the same optimal costs", {
  set.seed(1)
  y <- c(rnorm(100), rnorm(100, 2, 3), rnorm(100, -1, 0.5))
  pelt <- dust.meanVar(y, method = "PELT", backend = "scalar")
  for (method in c("1D", "2D")) for (backend in c("scalar", "highway"))
  {
    res <- dust.meanVar(y, method = method, backend = backend)
    expect_equal(res$costQ, pelt$costQ, info = paste(method, backend))
    expect_equal(as.numeric(res$changepoints), as.numeric(pelt$changepoints))
  }
  expect_equal(pelt$changepoints, c(100, 200, 300), tolerance = 0.03)
})

test_that("dust.object.meanVar gives the same result with data added step by step", {
  set.seed(2)
  y <- c(rnorm(80), rnorm(80, 1, 2))
  for (backend in c("scalar", "highway"))
  {
    one <- dust.meanVar(y, 4 * log(160), "2D", backend)
    ob <- dust.object.meanVar("2D", backend)
    ob$append_data(y[1:50], 4 * log(160))
    ob$update_partition()
    ob$append_data(y[51:160], NULL)
    ob$update_partition()
    expect_equal(ob$get_partition(), one, info = backend)
  }
})

test_that("meanVar reproduces the original dust pruning sets", {
  y <- sin((1:40) * 1.31) + cos((1:40) * 0.37)
  # t = 4 to 6: the variance of a one-point segment is now exactly 0
  # (sums with rounding errors), an index is pruned one step sooner
  expected <- list(
    `1D` = c(1,2,3,3,4,4,5,5,5,6,7,8,9,9,9,7,7,6,6,7,7,8,9,9,9,9,8,7,8,9,8,7,7,8,8,8,8,8,9,10),
    `2D` = c(1,2,3,3,4,4,5,5,5,6,6,7,8,8,8,7,6,6,6,7,7,8,8,7,6,7,7,7,8,9,7,6,7,8,6,7,7,8,9,10))
  for (method in c("1D", "2D")) for (backend in c("scalar", "highway"))
  {
    res <- dust.meanVar(y, 4 * log(40), method, backend)
    expect_identical(as.integer(res$nb), as.integer(expected[[method]]))
    expect_identical(as.integer(res$lastIndexSet), c(40L, 39L, 38L, 37L, 36L, 35L, 33L, 32L, 29L, 28L, 0L))
  }
})

test_that("meanVar input errors", {
  expect_error(dust.meanVar(c(1, NA, 3)))
  expect_error(dust.meanVar(rnorm(10), method = "3D"))
  expect_error(dust.meanVar(rnorm(10), penalty = -1))
  ob <- dust.object.meanVar()
  ob$append_data(rnorm(10), 2)
  expect_error(ob$append_data(rnorm(10), 3))
})
