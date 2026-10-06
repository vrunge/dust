test_that("multivariate generation agrees with successive 1D calls for every model", {
  types <- c("gauss", "exp", "poisson", "geom", "bern", "binom", "negbin", "variance")
  parameters <- cbind(c(0.3, 0.7), c(0.8, 0.2))
  gamma <- cbind(c(1, 0.5), c(0.8, 1))
  for(type in types) {
    set.seed(42)
    actual <- dataGenerator_MD(c(3, 7), parameters, sdNoise = c(0.2, 0.6),
                               gamma = gamma, nbTrials = c(3, 8),
                               nbSuccess = c(2, 5), type = type)
    set.seed(42)
    expected <- rbind(
      dataGenerator_1D(c(3, 7), parameters[, 1], sdNoise = 0.2,
                       gamma = gamma[, 1], nbTrials = 3, nbSuccess = 2, type = type),
      dataGenerator_1D(c(3, 7), parameters[, 2], sdNoise = 0.6,
                       gamma = gamma[, 2], nbTrials = 8, nbSuccess = 5, type = type))
    expect_equal(actual, expected, info = type)
    expect_equal(dim(dataGenerator_MD(type = type)), c(2L, 100L), info = type)
  }
})

test_that("Gaussian decay broadcasts across components and restarts at each segment", {
  parameters <- data.frame(first = c(2, 4), second = c(6, 8))
  expected <- rbind(c(2, 2, 2, 4, 2, 1), c(6, 6, 6, 8, 4, 2))
  expect_equal(dataGenerator_MD(c(3, 6), parameters, sdNoise = 0,
                                gamma = c(1, 0.5)), expected)
  gamma <- data.frame(first = c(1, 0.5), second = c(0.5, 1))
  expected <- rbind(c(2, 2, 2, 4, 2, 1), c(6, 3, 1.5, 8, 8, 8))
  expect_equal(dataGenerator_MD(c(3, 6), parameters, sdNoise = 0, gamma = gamma),
               expected)
  expect_equal(dataGenerator_MD(c(3, 6), parameters, sdNoise = 0,
                                gamma = as.matrix(gamma)), expected)
  expect_equal(dataGenerator_MD(c(3, 6), parameters, sdNoise = 0,
                                gamma = matrix(c(1, 0.5), nrow = 1)),
               rbind(c(2, 2, 2, 4, 4, 4), c(6, 3, 1.5, 8, 4, 2)))
})

test_that("single components and single observations retain matrix dimensions", {
  expect_equal(dataGenerator_MD(c(2, 4), matrix(c(2, 5), ncol = 1), sdNoise = 0),
               matrix(c(2, 2, 5, 5), nrow = 1))
  expect_equal(dataGenerator_MD(1, matrix(c(2, 5), nrow = 1), sdNoise = 0),
               matrix(c(2, 5), ncol = 1))
  expect_equal(dataGenerator_MD(1, matrix(2), sdNoise = 0), matrix(2))
})

test_that("generated Gaussian matrices can be passed directly to dust.MD", {
  y <- dataGenerator_MD(c(6, 12), cbind(c(0, 5), c(0, -3)), sdNoise = 0)
  fit <- dust.MD(y, penalty = 1, model = "gauss", method = "exact", backend = "scalar")
  expect_equal(fit$changepoints, c(6, 12))
})

test_that("multivariate generation rejects incompatible shapes and invalid values", {
  parameters <- matrix(0.5, nrow = 2, ncol = 2)
  for(chpts in list(numeric(), c(2, 2), c(2, 1), c(0, 2), c(1.5, 3),
                   c(2, NA_real_), c(2, Inf), c("2", "4")))
    expect_error(dataGenerator_MD(chpts, parameters), "chpts")
  for(value in list(c(1, 2), matrix(0.5, 1, 2), matrix(numeric(), 2, 0),
                   matrix(NA_real_, 2, 2), matrix(Inf, 2, 2),
                   data.frame(a = c("a", "b"), b = c(1, 2))))
    expect_error(dataGenerator_MD(c(2, 4), value), "parameters")
  for(value in list(NA_character_, character(), c("gauss", "exp"), "typo"))
    expect_error(dataGenerator_MD(c(2, 4), parameters, type = value), "type")
  expect_error(dataGenerator_MD(c(2, 4), parameters, sdNoise = 1:3), "sdNoise")
  expect_error(dataGenerator_MD(c(2, 4), parameters, sdNoise = c(1, -1)), "sdNoise")
  expect_error(dataGenerator_MD(c(2, 4), parameters, gamma = c(1, 0)), "gamma")
  for(value in list(numeric(), c(1, NA_real_), rep(1, 3), matrix(1, 2, 1),
                   matrix(1, 3, 2), data.frame(a = c("a", "b"), b = c(1, 1))))
    expect_error(dataGenerator_MD(c(2, 4), parameters, gamma = value), "gamma")
  expect_error(dataGenerator_MD(c(2, 4), parameters, nbTrials = 1:3,
                                type = "binom"), "nbTrials")
  expect_error(dataGenerator_MD(c(2, 4), parameters, nbTrials = 1.5,
                                type = "binom"), "nbTrials")
  expect_error(dataGenerator_MD(c(2, 4), parameters, nbSuccess = 1:3,
                                type = "negbin"), "nbSuccess")
  expect_error(dataGenerator_MD(c(2, 4), parameters, nbSuccess = 0,
                                type = "negbin"), "nbSuccess")
  expect_error(dataGenerator_MD(2, matrix(0, 1, 2), type = "variance"),
               "standard deviations must be positive")
  expect_error(dataGenerator_MD(2, matrix(0, 1, 2), type = "exp"),
               "rates must be positive")
  expect_error(dataGenerator_MD(2, matrix(1.5, 1, 2), type = "bern"),
               "probabilities")
})
