test_that("dataGenerator_1D gives data in the model domain", {
  set.seed(1)
  expect_length(dataGenerator_1D(chpts = c(50, 120), parameters = c(0, 1)), 120)
  expect_true(all(dataGenerator_1D(c(50, 100), c(2, 5), type = "poisson") %% 1 == 0))
  expect_true(all(dataGenerator_1D(c(50, 100), c(1, 2), type = "exp") > 0))
  expect_true(all(dataGenerator_1D(c(50, 100), c(0.3, 0.6), type = "geom") >= 1))
  expect_true(all(dataGenerator_1D(c(50, 100), c(0.3, 0.6), type = "bern") %in% 0:1))
  expect_true(all(dataGenerator_1D(c(50, 100), c(0.3, 0.6), nbTrials = 5, type = "binom") <= 5))
  expect_error(dataGenerator_1D(c(50, 40), c(0, 1)))
  expect_error(dataGenerator_1D(c(50, 100), c(0, 1, 2)))
  expect_error(dataGenerator_1D(c(50, 100), c(0.3, 1.2), type = "bern"))
})

test_that("dataGenerator_MD is a matrix of 1D series", {
  parameters <- cbind(c(2, 5), c(4, 1), c(3, 3))
  set.seed(2)
  y <- dataGenerator_MD(chpts = c(30, 60), parameters = parameters, type = "poisson")
  set.seed(2)
  z <- rbind(dataGenerator_1D(c(30, 60), parameters[, 1], type = "poisson"),
             dataGenerator_1D(c(30, 60), parameters[, 2], type = "poisson"),
             dataGenerator_1D(c(30, 60), parameters[, 3], type = "poisson"))
  expect_equal(unname(y), unname(z))
  expect_equal(dim(dataGenerator_MD(chpts = 10, parameters = matrix(0, 1, 1))), c(1, 10))
})

test_that("data_normalization_1D", {
  set.seed(3)
  y <- dataGenerator_1D(c(1000, 2000), c(0, 1), sdNoise = 3)
  expect_equal(sdDiff(data_normalization_1D(y)), 1)
  expect_equal(mean(data_normalization_1D(rpois(100, 4) + 1, "poisson")), 1)
  counts <- rbinom(50, 7, 0.4)
  expect_equal(data_normalization_1D(counts, "binom", size = 7), counts / 7)
  expect_error(data_normalization_1D(counts, "binom"))
})

test_that("segmentation_Cost_1D is the cost found by dust.1D", {
  set.seed(4)
  for (model in c("gauss", "poisson", "variance"))
  {
    y <- dataGenerator_1D(c(100, 200), list(gauss = c(0, 2), poisson = c(2, 6), variance = c(1, 3))[[model]], type = model)
    res <- dust.1D(y, 2 * log(200), model = model)
    cost <- segmentation_Cost_1D(y, res$changepoints, model)
    if (model == "gauss") cost <- cost - sum(y^2)
    k <- length(res$changepoints) - 1
    expect_equal(tail(res$costQ, 1), cost + k * 2 * log(200), info = model)
  }
})
