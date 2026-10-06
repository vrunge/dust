md_search_models <- c("gauss", "poisson", "exp", "geom", "bern", "binom",
                      "negbin", "variance")

md_search_data <- function(model, dimension = 3L, n = 36L) {
  set.seed(2309)
  mu <- rep(c(0, 1.5, -0.5), length.out = n)
  prob <- rep(c(0.2, 0.7, 0.4), each = n / 3)
  rows <- lapply(seq_len(dimension), function(i) switch(model,
    gauss = rnorm(n, mu), poisson = rpois(n, 2 + 4 * prob),
    exp = rexp(n, prob), geom = rgeom(n, prob) + 1,
    bern = rbinom(n, 1, prob), binom = rbinom(n, 10, prob) / 10,
    negbin = rnbinom(n, size = 10, prob = prob) / 10,
    variance = rnorm(n, sd = 3 * prob)))
  do.call(rbind, rows)
}

test_that("new MD searches match every PELT prefix cost and its changepoints", {
  for (model in md_search_models) for (dimension in c(1L, 3L)) {
    data <- md_search_data(model, dimension)
    for (penalty in c(0, 2, 2 * dimension * log(ncol(data)))) {
      ref <- dust.MD(data, penalty, model, "PELT", "scalar")
      op <- dust.MD(data, penalty, model, "OP", "scalar")
      expect_equal(ref$costQ, op$costQ, tolerance = 1e-9, info = model)
      expect_identical(ref$changepoints, op$changepoints, info = model)
      for (p in unique(c(1L, dimension))) for (loops in c(1L, 10L, 40L)) {
        for (method in c("iterative", "QN")) {
          info <- paste(model, dimension, penalty, p, loops, method)
          fit <- dust.MD(data, penalty, model, method, "scalar", p, loops)
          expect_equal(fit$costQ, ref$costQ, tolerance = 1e-9, info = info)
          expect_identical(fit$changepoints, ref$changepoints, info = info)
          expect_true(all(fit$nb <= ref$nb), info = info)
        }
      }
    }
  }
})

test_that("new MD searches handle boundary faces, ties and dependent constraints", {
  for (model in md_search_models) {
    value <- switch(model, gauss = 0, poisson = 0, exp = 1, geom = 1,
                     bern = 0, binom = 1, negbin = 0, variance = 1)
    data <- md_search_data(model)
    # A fixed boundary row and two identical variable rows give singular
    # curvature and exercise optimization on a lower dimensional face.
    cases <- list(matrix(value, 3, 18), rbind(value, data[1, ], data[1, ]))
    for (data in cases) for (penalty in c(0, 1)) {
      ref <- dust.MD(data, penalty, model, "PELT", "scalar")
      for (method in c("iterative", "QN")) {
        fit <- dust.MD(data, penalty, model, method, "scalar", constraints = 2)
        expect_equal(fit$costQ, ref$costQ, tolerance = 1e-9,
                     info = paste(model, method, penalty))
        expect_identical(fit$changepoints, ref$changepoints,
                         info = paste(model, method, penalty))
        expect_true(all(fit$nb <= ref$nb))
      }
    }
  }
})

test_that("new MD methods resume identically for all models and backends", {
  for (model in md_search_models) for (method in c("iterative", "QN")) {
    data <- md_search_data(model, dimension = 6L)
    for (backend in c("scalar", "highway")) {
      one <- dust.MD(data, 3, model, method, backend, constraints = 3)
      ref <- dust.MD(data, 3, model, "PELT", backend)
      ob <- dust.object.MD(model, method, backend, constraints = 3)
      ob$append_data(data[, 1:13, drop = FALSE], 3)
      ob$update_partition()
      ob$append_data(data[, 14:36, drop = FALSE], NULL)
      ob$update_partition()
      fit <- ob$get_partition()
      expect_equal(fit, one, tolerance = 1e-10, info = paste(model, method, backend))
      expect_equal(fit$costQ, ref$costQ, tolerance = 1e-9)
      expect_identical(fit$changepoints, ref$changepoints)
      expect_identical(ob$get_info()$pruning_algo, method)
    }
  }
})

test_that("both searches find additional pruning witnesses for every model", {
  for (model in md_search_models) {
    data <- md_search_data(model)
    penalty <- 2 * nrow(data) * log(ncol(data))
    ref <- dust.MD(data, penalty, model, "PELT", "scalar")
    for (method in c("iterative", "QN")) {
      fit <- dust.MD(data, penalty, model, method, "scalar", 3, 40)
      expect_lt(sum(fit$nb), sum(ref$nb), label = paste(model, method))
      expect_equal(fit$costQ, ref$costQ, tolerance = 1e-9)
      expect_identical(fit$changepoints, ref$changepoints)
    }
  }
})

test_that("epsilon uses a 1000-iteration cap and explicit nbIterations wins", {
  for (method in c("coordinateDescent", "iterative", "QN")) {
    epsilon_object <- dust.object.MD(method = method, epsilon = 1e9)
    expect_identical(epsilon_object$get_info()$nbIterations, 1000L)
    expect_identical(epsilon_object$get_info()$epsilon, 1e9)
    fixed_object <- dust.object.MD(method = method, nbIterations = 3L,
                                   epsilon = 1e9)
    expect_identical(fixed_object$get_info()$nbIterations, 3L)
    expect_null(fixed_object$get_info()$epsilon)
  }
  expect_identical(dust.object.MD()$get_info()$nbIterations, 10L)
  expect_null(dust.object.MD()$get_info()$epsilon)
  expect_null(dust.object.MD(method = "exact", epsilon = 1e-6)$get_info()$epsilon)
  expect_error(dust.object.MD(epsilon = -1), "epsilon")
  expect_error(dust.object.MD(epsilon = Inf), "epsilon")
  expect_error(dust.object.MD(method = "randomEval", epsilon = 1e-6), "epsilon")
})

test_that("epsilon stopping preserves valid pruning across MD models", {
  for (model in md_search_models) {
    data <- md_search_data(model)
    penalty <- 2 * nrow(data) * log(ncol(data))
    reference <- dust.MD(data, penalty, model, "PELT", "scalar")
    for (method in c("coordinateDescent", "iterative", "QN")) {
      epsilon_fit <- dust.MD(data, penalty, model, method, "scalar",
                             epsilon = 1e9)
      fixed_fit <- dust.MD(data, penalty, model, method, "scalar",
                           nbIterations = 1L)
      priority_fit <- dust.MD(data, penalty, model, method, "scalar",
                              nbIterations = 1L, epsilon = 1e9)
      expect_identical(priority_fit, fixed_fit,
                       info = paste(model, method))
      expect_identical(epsilon_fit$nb, fixed_fit$nb,
                       info = paste(model, method))
      expect_equal(epsilon_fit$costQ, reference$costQ, tolerance = 1e-9,
                   info = paste(model, method))
      expect_identical(epsilon_fit$changepoints, reference$changepoints,
                       info = paste(model, method))
      expect_true(all(epsilon_fit$nb <= reference$nb),
                  info = paste(model, method))
    }
  }
})
