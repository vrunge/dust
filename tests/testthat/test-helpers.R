test_that("the exported cost helpers use the model's statistic and validate endpoints", {
  y <- c(-2, 1, -3, 2)
  expect_equal(segmentation_Cost_1D(y, 4, "variance"),
               length(y) * (1 + log(mean(y^2))) / 2)
  expect_error(Cost_1D(c(0, 1, 2), 1, 3, "typo"), "unknown model")
  expect_error(Cost_1D(c(0, 1, 2), 2, 2, "gauss"), "nonempty segment")
  expect_error(segmentation_Cost_1D(y, integer(), "gauss"), "chpts")
  expect_error(segmentation_Cost_1D(y, c(4, 2), "gauss"), "chpts")
  expect_error(segmentation_Cost_1D(c(0, 2, 3), 3, "binom"), "domain")
})

test_that("normalization uses known count-model sizes and handles zero data", {
  expect_equal(data_normalization_1D(c(0, 2, 3), "binom", size = 5),
               c(0, 0.4, 0.6))
  expect_equal(data_normalization_1D(c(0, 2, 3), "negbin", size = 10),
               c(0, 0.2, 0.3))
  expect_error(data_normalization_1D(c(0, 2), "binom"), "size")
  expect_error(data_normalization_1D(c(0, 2), "negbin"), "size")
  expect_equal(data_normalization_1D(rep(0, 5), "poisson"), rep(0, 5))
  expect_error(data_normalization_1D(rep(0, 5), "gauss"), "positive")
  expect_error(data_normalization_1D(c(0, 1), "exp"), "strictly positive")
  expect_error(data_normalization_1D(c(1, 2, 3), "variance"), "zero residual")
  expect_error(data_normalization_1D(c(0, 1), "geom"), "at least one")
  expect_error(data_normalization_1D(c(0, 2), "bern"), "Bernoulli")
  expect_error(sdDiff(1:5, "typo"), "method")
  expect_equal(sdDiff(1:5), sdDiff(1:5, "HALL"))
})

test_that("the Gaussian generator accepts mixed decay factors", {
  expect_equal(dataGenerator_1D(c(3, 6), c(2, 4), sdNoise = 0,
                                gamma = c(1, 0.5), type = "gauss"),
               c(2, 2, 2, 4, 2, 1))
})
