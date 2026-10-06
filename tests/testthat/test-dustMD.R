md_fixture <- function(model) {
  levels <- switch(model,
    gauss = c(0, 3, -1), poisson = c(1, 6, 2),
    exp = c(1, 0.25, 2), geom = c(2, 6, 3),
    bern = c(0.1, 0.9, 0.2), binom = c(0.1, 0.9, 0.2),
    negbin = c(0.5, 4, 1), variance = c(1, 4, 2))
  first <- rep(levels, each = 6)
  rbind(first, rev(first), first * if (model %in% c("bern", "binom")) 0.9 else 1.1)
}

test_that("MD pruning methods preserve the OP optimum for all models", {
  methods <- c("PELT", "randomEval", "coordinateDescent", "iterative", "QN")
  models <- c("gauss", "poisson", "exp", "geom", "bern", "binom",
              "negbin", "variance")
  for (model in models) {
    y <- md_fixture(model)
    penalty <- 2 * nrow(y) * log(ncol(y))
    reference <- dust.MD(y, penalty, model, "OP", "scalar")
    expect_identical(reference$backend, "scalar", info = model)
    expect_equal(reference$nb, seq_len(ncol(y)), info = model)
    expect_equal(reference$lastIndexSet, ncol(y):0, info = model)
    for (method in methods) {
      set.seed(123)
      fit <- dust.MD(y, penalty, model, method, "scalar", constraints = 3)
      expect_equal(fit$costQ, reference$costQ, tolerance = 1e-8,
                   info = paste(model, method))
      expect_equal(fit$changepoints, reference$changepoints,
                   info = paste(model, method))
      expect_equal(length(fit$nb), ncol(y), info = model)
      expect_lte(length(fit$lastIndexSet), length(reference$lastIndexSet))
    }
  }
})

test_that("MD retains the OP optimum on irregular deterministic series", {
  wave <- sin(seq_len(24) * 0.7) + cos(seq_len(24) * 1.3) / 4
  observations <- list(
    gauss = wave,
    poisson = round(pmax(0, 3 + 2 * wave)),
    exp = exp(wave / 2),
    geom = pmax(1, round(3 + wave)),
    bern = as.numeric(wave > 0),
    binom = pmin(1, pmax(0, 0.5 + wave / 4)),
    negbin = pmax(0, 2 + wave),
    variance = exp(wave / 2))
  for (model in names(observations)) {
    values <- observations[[model]]
    y <- rbind(values, rev(values), values)
    for (penalty in c(0, 2.5)) {
      reference <- dust.MD(y, penalty, model, "OP", "scalar")
      for (p in c(1, nrow(y))) {
        for (method in c("PELT", "coordinateDescent", "randomEval", "iterative", "QN")) {
          set.seed(44)
          fit <- dust.MD(y, penalty, model, method, "scalar", p)
          expect_equal(fit$costQ, reference$costQ, tolerance = 1e-8,
                       info = paste(model, penalty, p, method))
          expect_equal(fit$changepoints, reference$changepoints,
                       info = paste(model, penalty, p, method))
        }
      }
    }
  }
})

test_that("MD uses available nearest predecessors before the constraint limit", {
  y <- md_fixture("gauss")
  pen <- 1
  pelt <- dust.MD(y, pen, method = "PELT", backend = "scalar")
  for (p in seq_len(nrow(y))) {
    fit <- dust.MD(y, pen, method = "coordinateDescent",
                   backend = "scalar", constraints = p)
    expect_equal(fit$costQ, pelt$costQ)
    expect_identical(fit$changepoints, pelt$changepoints)
  }
  expect_equal(dust.object.MD()$get_info()$constraints, 0)

  set.seed(20261006)
  y <- matrix(rnorm(200), nrow = 10)
  penalty <- 20 * log(1000)
  exact <- dust.MD(y, penalty, method = "exact", backend = "scalar",
                   constraints = 10)
  expect_equal(exact$nb[2], 1)
  for (method in c("coordinateDescent", "iterative", "QN")) {
    fit <- dust.MD(y, penalty, method = method, backend = "scalar",
                   constraints = 10, nbIterations = 100)
    expect_equal(fit$nb[2], exact$nb[2], info = method)
    expect_equal(fit$costQ, exact$costQ, tolerance = 1e-9, info = method)
    expect_identical(fit$changepoints, exact$changepoints, info = method)
  }
})

test_that("MD object resumes from matrix batches", {
  y <- md_fixture("gauss")
  pen <- 2 * nrow(y) * log(ncol(y))
  for (method in c("OP", "PELT", "coordinateDescent", "randomEval", "iterative", "QN")) {
    set.seed(91)
    one <- dust.MD(y, pen, method = method, backend = "scalar", constraints = 2)
    set.seed(91)
    ob <- dust.object.MD(method = method, backend = "scalar", constraints = 2)
    ob$append_data(y[, 1:5, drop = FALSE], pen)
    expect_error(ob$get_partition(), "update_partition")
    ob$update_partition()
    ob$append_data(y[, 6:11, drop = FALSE], NULL)
    ob$update_partition()
    ob$append_data(y[, 12:18, drop = FALSE], NULL)
    ob$update_partition()
    fit <- ob$get_partition()
    expect_equal(fit$changepoints, one$changepoints, info = method)
    expect_equal(fit$costQ, one$costQ, tolerance = 1e-9, info = method)
    expect_equal(ob$get_info()$constraints, 2)
    expect_equal(ob$get_info()$dimension, nrow(y))
  }
})

test_that("MD default penalty scales with dimension", {
  y <- c(rep(0, 5), rep(3, 5), rep(-1, 5))
  multi <- rbind(y, y, y)
  single <- dust.1D(y, method = "OP", backend = "scalar")
  fit <- dust.MD(multi, method = "OP", backend = "scalar")
  expect_equal(fit$costQ, 3 * single$costQ, tolerance = 1e-9)
  expect_equal(fit$changepoints, single$changepoints)
})

test_that("one-row MD OP has the 1D objective for every model", {
  for (model in c("gauss", "poisson", "exp", "geom", "bern", "binom",
                  "negbin", "variance")) {
    y <- md_fixture(model)[1, ]
    pen <- 2 * log(length(y))
    one <- dust.1D(y, pen, model, "OP", "scalar")
    multi <- dust.MD(matrix(y, nrow = 1), pen, model, "OP", "scalar")
    expect_equal(multi$costQ, one$costQ, tolerance = 1e-8, info = model)
    expect_equal(multi$changepoints, one$changepoints, info = model)
  }
})

test_that("MD handles model boundaries and validates its inputs", {
  boundary <- list(
    poisson = rbind(c(0, 0, 2, 2), c(1, 1, 0, 0)),
    geom = rbind(c(1, 1, 4, 4), c(3, 3, 1, 1)),
    bern = rbind(c(0, 0, 1, 1), c(1, 1, 0, 0)),
    binom = rbind(c(0, 0, 1, 1), c(1, 1, 0, 0)),
    negbin = rbind(c(0, 0, 3, 3), c(2, 2, 0, 0)))
  for (model in names(boundary)) {
    ref <- dust.MD(boundary[[model]], 1, model, "OP", "scalar")
    for (method in c("coordinateDescent", "iterative", "QN")) {
      fit <- dust.MD(boundary[[model]], 1, model, method, "scalar", constraints = 2)
      expect_equal(fit$costQ, ref$costQ, info = paste(model, method))
      expect_equal(fit$changepoints, ref$changepoints, info = paste(model, method))
    }
  }
  expect_error(dust.MD(1:3), "numeric matrix")
  expect_error(dust.MD(matrix(1:4, 2), constraints = 3), "constraints")
  expect_error(dust.MD(matrix(c(1, NA), 1)), "model domain")
  expect_error(dust.MD(matrix(c(0, 1), 1), model = "exp"), "model domain")
  expect_error(dust.object.MD(nbIterations = 0), "nbIterations")
  expect_error(dust.object.MD(method = "DUST"))
  ob <- dust.object.MD()
  ob$append_data(matrix(1:4, 2), 1)
  expect_error(ob$append_data(matrix(1:3, 1), 1), "number of rows")
  expect_error(ob$append_data(matrix(1:2, 2), 2), "penalty cannot change")
})

test_that("MD random evaluation is reproducible under set.seed", {
  y <- md_fixture("gauss")
  set.seed(765)
  first <- dust.MD(y, penalty = 1, method = "randomEval", backend = "scalar")
  set.seed(765)
  second <- dust.MD(y, penalty = 1, method = "randomEval", backend = "scalar")
  expect_identical(first, second)
})

test_that("Gaussian exact joint pruning preserves the OP solution", {
  set.seed(374)
  series <- list(
    md_fixture("gauss"),
    rbind(rnorm(48), rnorm(48), rnorm(48)),
    rbind(rep(0, 24), rep(0, 24)),
    rbind(seq_len(24) / 10, rep(1, 24)))
  for (y in series) {
    for (penalty in c(0, 1, 2 * nrow(y) * log(ncol(y)))) {
      reference <- dust.MD(y, penalty, method = "OP", backend = "scalar")
      pelt <- dust.MD(y, penalty, method = "PELT", backend = "scalar")
      for (constraints in seq_len(nrow(y))) {
        fit <- dust.MD(y, penalty, method = "exact", backend = "scalar",
                       constraints = constraints)
        expect_equal(fit$costQ, reference$costQ, tolerance = 1e-9)
        expect_equal(fit$changepoints, reference$changepoints)
        expect_true(all(fit$nb <= pelt$nb))
      }
    }
  }
})

test_that("Gaussian exact works across appends and other models use PELT", {
  y <- md_fixture("gauss")
  one <- dust.MD(y, 3, method = "exact", backend = "scalar", constraints = 2)
  ob <- dust.object.MD(method = "exact", backend = "scalar", constraints = 2)
  ob$append_data(y[, 1:7, drop = FALSE], 3)
  ob$update_partition()
  ob$append_data(y[, 8:18, drop = FALSE], NULL)
  ob$update_partition()
  expect_equal(ob$get_partition(), one)
  expect_identical(ob$get_info()$pruning_algo, "exact")
  for (model in c("poisson", "exp", "geom", "bern", "binom", "negbin", "variance")) {
    y <- md_fixture(model)
    pelt <- dust.MD(y, 2, model, "PELT", "scalar", constraints = 2)
    exact <- dust.MD(y, 2, model, "exact", "scalar", constraints = 2)
    expect_identical(exact, pelt, info = model)
    object <- dust.object.MD(model, "exact", "scalar", constraints = 2)
    expect_identical(object$get_info()$pruning_algo, "PELT", info = model)
    object$append_data(y[, 1:7, drop = FALSE], 2)
    object$update_partition()
    object$append_data(y[, 8:18, drop = FALSE], NULL)
    object$update_partition()
    expect_identical(object$get_partition(), pelt, info = model)
  }
})

test_that("Gaussian exact method improves pruning and agrees across backends", {
  set.seed(123)
  y <- rbind(rnorm(100), rnorm(100))
  pelt <- dust.MD(y, 4, method = "PELT", backend = "scalar")
  scalar <- dust.MD(y, 4, method = "exact", backend = "scalar")
  vector <- dust.MD(y, 4, method = "exact", backend = "highway")
  expect_lt(sum(scalar$nb), sum(pelt$nb))
  expect_equal(vector$costQ, scalar$costQ, tolerance = 1e-9)
  expect_equal(vector$changepoints, scalar$changepoints)
  expect_identical(vector$backend, DUST.1D.HW.backend())
})

test_that("joint Gaussian exact uses multiple constraints to prune further", {
  set.seed(1)
  y <- rbind(rnorm(80), rnorm(80))
  one <- dust.MD(y, 5, method = "exact", backend = "scalar", constraints = 1)
  joint <- dust.MD(y, 5, method = "exact", backend = "scalar", constraints = 2)
  reference <- dust.MD(y, 5, method = "PELT", backend = "scalar")
  expect_lt(sum(joint$nb), sum(one$nb))
  expect_equal(joint$costQ, reference$costQ, tolerance = 1e-9)
  expect_identical(joint$changepoints, reference$changepoints)
})

test_that("MD Highway and scalar return the same optimum", {
  for (model in c("gauss", "poisson", "exp", "geom", "bern", "binom",
                  "negbin", "variance")) {
    y <- rbind(md_fixture(model), md_fixture(model))
    for (method in c("OP", "PELT", "coordinateDescent", "randomEval", "iterative", "QN", "exact")) {
      set.seed(21)
      scalar <- dust.MD(y, penalty = 8, model = model, method = method,
                        backend = "scalar", constraints = 3)
      set.seed(21)
      vector <- dust.MD(y, penalty = 8, model = model, method = method,
                        backend = "highway", constraints = 3)
      expect_equal(vector$costQ, scalar$costQ, tolerance = 1e-8,
                   info = paste(model, method))
      expect_equal(vector$changepoints, scalar$changepoints,
                   info = paste(model, method))
      expect_identical(vector$backend, DUST.1D.HW.backend())
    }
  }
})

test_that("MD Highway candidate arrays resume across uneven batches", {
  set.seed(700)
  n <- 360L
  level <- rep(c(1, 4, 2), each = n / 3L)
  series <- list(
    gauss = matrix(rnorm(3L * n, rep(level, each = 3L)), 3L),
    poisson = matrix(rpois(3L * n, rep(level, each = 3L)), 3L)
  )
  for (model in names(series)) {
    y <- series[[model]]
    for (method in c("OP", "PELT", "coordinateDescent", "iterative", "QN", "exact")) {
      one <- dust.MD(y, penalty = 10, model = model, method = method,
                     backend = "highway", constraints = 3L)
      object <- dust.object.MD(model = model, method = method,
                               backend = "highway", constraints = 3L)
      for (indices in list(1:41, 42:137, 138:360)) {
        object$append_data(y[, indices, drop = FALSE], 10)
        object$update_partition()
      }
      expect_identical(object$get_partition(), one,
                       info = paste(model, method))
    }
  }
})

test_that("binomial MD keeps a rounded endpoint in the segment domain", {
  y <- matrix(c(.2, .5, .4, .6, .7, .8, .4, .1, .2,
                .3, .2, 0, .8, .7, 1), nrow = 3)
  reference <- dust.MD(y, 0, "binom", "OP", "highway")
  for (method in c("OP", "PELT", "coordinateDescent", "randomEval",
                   "iterative", "QN", "exact"))
    for (backend in c("scalar", "highway")) {
      fit <- dust.MD(y, 0, "binom", method, backend, constraints = 3)
      expect_equal(fit$costQ, reference$costQ, tolerance = 1e-12,
                   info = paste(method, backend))
    }
})
