### 2-row data with 2 changes for each model
data_MD <- function(model, n = 150)
{
  parameters <- list(gauss = c(0, 2, -1), poisson = c(2, 6, 3), exp = c(1, 0.3, 2),
                     geom = c(0.6, 0.2, 0.5), bern = c(0.2, 0.8, 0.4), binom = c(0.2, 0.7, 0.4),
                     negbin = c(0.6, 0.2, 0.4), variance = c(1, 3, 0.5))[[model]]
  y <- dataGenerator_MD(chpts = round(n * (1:3) / 3), parameters = cbind(parameters, rev(parameters)),
                        nbTrials = 10, nbSuccess = 10, type = model)
  if (model %in% c("binom", "negbin")) y <- y / 10
  y
}
models <- c("gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance")


test_that("every dust.MD method gives the optimal segmentation (as OP)", {
  set.seed(1)
  for (model in models)
  {
    y <- data_MD(model)
    op <- dust.MD(y, model = model, method = "OP")
    for (method in c("exact", "coordinateDescent", "QN", "randomEval", "PELT"))
      for (constraints in 1:2)
      {
        res <- dust.MD(y, model = model, method = method, constraints = constraints)
        info <- paste(model, method, constraints)
        expect_equal(res$costQ, op$costQ, tolerance = 1e-9, info = info)
        expect_equal(res$changepoints, op$changepoints, info = info)
      }
  }
})

test_that("dust.MD with one row is dust.1D", {
  set.seed(2)
  y <- dataGenerator_1D(chpts = c(100, 200), parameters = c(2, 5), type = "poisson")
  md <- dust.MD(matrix(y, nrow = 1), 2 * log(200), model = "poisson", method = "OP")
  expect_equal(md$costQ, dust.1D(y, 2 * log(200), model = "poisson", method = "OP")$costQ)
})

test_that("Gaussian exact prunes as the maximum of the decision function (converged QN)", {
  set.seed(3)
  y <- matrix(rnorm(3 * 400), nrow = 3)
  for (constraints in 1:3)
  {
    exact <- dust.MD(y, method = "exact", constraints = constraints)
    qn <- dust.MD(y, method = "QN", constraints = constraints, nbIterations = 1000)
    expect_identical(exact$nb, qn$nb, info = constraints)
  }
  expect_lt(sum(dust.MD(y, method = "exact", constraints = 3)$nb), sum(dust.MD(y, method = "PELT")$nb) / 3)
})

test_that("dust.object.MD gives the same result with data added step by step", {
  set.seed(4)
  y <- data_MD("gauss")
  for (method in c("exact", "coordinateDescent", "QN"))
  {
    one <- dust.MD(y, 3 * log(150), method = method, constraints = 2)
    ob <- dust.object.MD(method = method, constraints = 2)
    ob$append_data(y[, 1:40, drop = FALSE], 3 * log(150))
    ob$update_partition()
    ob$append_data(y[, 41:150], NULL)
    ob$update_partition()
    expect_equal(ob$get_partition(), one, info = method)
  }
})

test_that("dust.MD input errors", {
  y <- matrix(rnorm(20), nrow = 2)
  expect_error(dust.MD(rnorm(10)))
  expect_error(dust.MD(y, constraints = 3))
  expect_error(dust.MD(y, method = "iterative"))
  expect_error(dust.MD(y, nbIterations = 0))
  expect_error(dust.MD(y, method = "randomEval", epsilon = 1e-6))
  expect_error(dust.MD(y, model = "poisson"))
  ob <- dust.object.MD()
  ob$append_data(y, 1)
  expect_error(ob$append_data(matrix(1, 3, 2), 1))
})
