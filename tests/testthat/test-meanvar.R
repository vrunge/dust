test_that("meanVar methods give the same optimal costs", {
  set.seed(1)
  y <- c(rnorm(100), rnorm(100, 2, 3), rnorm(100, -1, 0.5))
  pelt <- dust.meanVar(y, method = "PELT")
  for (method in c("1D", "2D"))
  {
    res <- dust.meanVar(y, method = method)
    expect_equal(res$costQ, pelt$costQ, info = method)
    expect_equal(as.numeric(res$changepoints), as.numeric(pelt$changepoints))
  }
  expect_equal(pelt$changepoints, c(100, 200, 300), tolerance = 0.03)
})

test_that("dust.object.meanVar gives the same result with data added step by step", {
  set.seed(2)
  y <- c(rnorm(80), rnorm(80, 1, 2))
 
  {
    one <- dust.meanVar(y, 3 * log(160), "2D")
    ob <- dust.object.meanVar("2D")
    ob$append_data(y[1:50], 3 * log(160))
    ob$update_partition()
    ob$append_data(y[51:160], NULL)
    ob$update_partition()
    expect_equal(ob$get_partition(), one)
  }
})

test_that("meanVar DUST keeps the optimal costs and fewer indices than PELT", {
  y <- sin((1:40) * 1.31) + cos((1:40) * 0.37)
  pelt <- dust.meanVar(y, 3 * log(40), "PELT")
  for (method in c("1D", "2D"))
  {
    res <- dust.meanVar(y, 3 * log(40), method)
    expect_equal(res$costQ, pelt$costQ, info = method)
    expect_true(all(res$nb <= pelt$nb), info = method)
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
