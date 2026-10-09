### data with 3 changes for each model (normalized)
data_1D <- function(model, n = 300)
{
  parameters <- list(gauss = c(0, 2, -1), poisson = c(2, 6, 3), exp = c(1, 0.3, 2),
                     geom = c(0.6, 0.2, 0.5), bern = c(0.2, 0.8, 0.4), binom = c(0.2, 0.7, 0.4),
                     negbin = c(0.6, 0.2, 0.4), variance = c(1, 3, 0.5))[[model]]
  chpts <- round(n * (1:3) / 3)
  y <- dataGenerator_1D(chpts = chpts, parameters = parameters, nbTrials = 10, nbSuccess = 10, type = model)
  if (model %in% c("binom", "negbin")) data_normalization_1D(y, model, size = 10) else data_normalization_1D(y, model)
}
models <- c("gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance")


test_that("dust.1D gives the optimal segmentation (as OP) for every model and method", {
  set.seed(1)
  for (model in models)
  {
    y <- data_1D(model)
    op <- dust.1D(y, model = model, method = "OP")
    for (method in c("DUST", "DUSTib", "PELT"))
    {
      res <- dust.1D(y, model = model, method = method)
      info <- paste(model, method)
      expect_equal(res$costQ, op$costQ, tolerance = 1e-10, info = info)
      expect_equal(as.numeric(res$changepoints), as.numeric(op$changepoints), info = info)
      expect_true(all(res$nb <= op$nb), info = info)
    }
  }
})

test_that("dust.1D finds the change points", {
  set.seed(2)
  for (model in c("gauss", "poisson", "variance"))
  {
    y <- data_1D(model, n = 1500)
    expect_equal(dust.1D(y, model = model)$changepoints, c(500, 1000, 1500), tolerance = 0.02, info = model)
  }
})

test_that("DUST keeps few indices when there is no change", {
  set.seed(3)
  y <- rnorm(3000)
  expect_equal(dust.1D(y)$changepoints, 3000)
  expect_lt(mean(dust.1D(y)$nb), mean(dust.1D(y, method = "PELT")$nb) / 20)
})

test_that("boundary data (zeros, ones) do not break the pruning tests", {
  set.seed(5)
  y <- list(poisson = c(rep(0, 50), rpois(50, 3), rep(0, 30)),
            bern = c(rep(0, 40), rep(1, 40), rbinom(60, 1, 0.5)),
            geom = c(rep(1, 50), rgeom(50, 0.3) + 1))
  for (model in names(y)) for (method in c("DUST", "DUSTib"))
  {
    op <- dust.1D(y[[model]], model = model, method = "OP")
    res <- dust.1D(y[[model]], model = model, method = method)
    expect_equal(res$costQ, op$costQ, tolerance = 1e-10, info = paste(model, method))
  }
})

test_that("dust.object.1D gives the same result with data added step by step", {
  set.seed(4)
  for (model in c("gauss", "negbin", "variance"))
  {
    y <- data_1D(model)
    one <- dust.1D(y, 2 * log(300), model = model)
    ob <- dust.object.1D(model = model)
    ob$append_data(y[1:70], 2 * log(300))
    ob$update_partition()
    ob$append_data(y[71:300], NULL)
    ob$update_partition()
    res <- ob$get_partition()
    expect_equal(res$costQ, one$costQ, info = model)
    expect_equal(as.numeric(res$changepoints), as.numeric(one$changepoints))
  }
})

test_that("dust.1D input errors", {
  expect_error(dust.1D(rnorm(10), model = "gaus"))
  expect_error(dust.1D(rnorm(10), method = "dust"))
  expect_error(dust.1D(c(1, NA)))
  expect_error(dust.1D(numeric(0)))
  expect_error(dust.1D(rnorm(10), penalty = -1))
  expect_error(dust.1D(c(1, -1), model = "poisson"))
  expect_error(dust.1D(c(1, 0), model = "exp"))
  expect_error(dust.object.1D()$get_partition())
})

test_that("PELTpar gives the optimal segmentation (as OP), serial and threaded", {
  set.seed(11)
  for (model in models)
  {
    y <- data_1D(model, n = 600)
    op <- dust.1D(y, model = model, method = "OP")
    res <- dust.1D(y, model = model, method = "PELTpar")
    expect_equal(res$costQ, op$costQ, tolerance = 1e-10, info = model)
    expect_equal(as.numeric(res$changepoints), as.numeric(op$changepoints), info = model)
  }
  y <- rnorm(30000)  # no change: large candidate range, threaded scan
  pelt <- dust.1D(y, method = "PELT")
  for (threads in c(1L, 4L))
  {
    res <- dust.1D(y, method = "PELTpar", threads = threads)
    expect_equal(res$costQ, pelt$costQ, tolerance = 1e-10)
    expect_equal(res$changepoints, pelt$changepoints)
  }
})
