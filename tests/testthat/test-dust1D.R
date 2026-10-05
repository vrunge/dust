test_that("dust.1D recovers known change points for every cost model", {
  set.seed(42)
  chpts <- c(100, 200, 300)

  cases <- list(
    gauss    = list(y = dataGenerator_1D(chpts, c(0, 3, 0), sdNoise = 1, type = "gauss"), normalize = TRUE),
    poisson  = list(y = dataGenerator_1D(chpts, c(2, 8, 2), type = "poisson"), normalize = FALSE),
    exp      = list(y = dataGenerator_1D(chpts, c(1, 4, 1), type = "exp"), normalize = FALSE),
    geom     = list(y = dataGenerator_1D(chpts, c(0.7, 0.2, 0.7), type = "geom"), normalize = FALSE),
    bern     = list(y = dataGenerator_1D(chpts, c(0.2, 0.8, 0.2), type = "bern"), normalize = FALSE),
    binom    = list(y = dataGenerator_1D(chpts, c(0.2, 0.8, 0.2), nbTrials = 10, type = "binom") / 10, normalize = FALSE),
    negbin   = list(y = dataGenerator_1D(chpts, c(0.7, 0.2, 0.7), nbSuccess = 10, type = "negbin") / 10, normalize = FALSE),
    variance = list(y = dataGenerator_1D(chpts, c(1, 5, 1), type = "variance"), normalize = FALSE)
  )

  # Fixed reference changepoints: a regression guard, not a statement that
  # these are "the" right answer -- what matters here is that the result
  # does not silently change when the C++ engine is refactored.
  reference <- list(
    gauss    = c(100, 200, 300),
    poisson  = c(100, 200, 300),
    exp      = c(96, 202, 300),
    geom     = c(100, 200, 300),
    bern     = c(101, 202, 300),
    binom    = c(100, 200, 300),
    negbin   = c(100, 200, 300),
    variance = c(100, 199, 300)
  )

  for (model in names(cases)) {
    y <- cases[[model]]$y
    if (cases[[model]]$normalize) y <- data_normalization_1D(y, type = model)
    res <- dust.1D(data = y, model = model, penalty = 2 * log(length(y)), backend = "scalar")
    expect_identical(res$changepoints, reference[[model]], info = model)
  }
})

test_that("all dualmax_algo variants run and agree on the segmentation", {
  set.seed(7)
  y <- dataGenerator_1D(c(100, 200, 300), c(0, 3, 0), sdNoise = 1, type = "gauss")
  y <- data_normalization_1D(y, type = "gauss")

  methods <- c("DUST", "DUSTib", "PELT", "OP")

  expected <- c(100, 201, 300)
  for (m in methods) {
    object <- dust.object.1D(model = "gauss", method = m, backend = "scalar")
    res <- object$dust(y, 2 * log(length(y)))
    expect_identical(unlist(as.list(res$changepoints)), expected, info = m)
  }
})

test_that("dust.1D OP engine never prunes and matches PELT's costQ exactly, for every cost model", {
  # Regression guard for the DUST_1D_OP_T<Model> engine (1D_OP_Impl.h):
  # the old method="OP" silently pruned a data-dependent amount (routed
  # through the general DualMaxPolicy engine), so this both checks
  # ground-truth agreement with PELT and that nb is trivially t at every
  # step (i.e. the active set is never reduced).
  set.seed(42)
  chpts <- c(100, 200, 300)

  cases <- list(
    gauss    = list(y = dataGenerator_1D(chpts, c(0, 3, 0), sdNoise = 1, type = "gauss"), normalize = TRUE),
    poisson  = list(y = dataGenerator_1D(chpts, c(2, 8, 2), type = "poisson"), normalize = FALSE),
    exp      = list(y = dataGenerator_1D(chpts, c(1, 4, 1), type = "exp"), normalize = FALSE),
    geom     = list(y = dataGenerator_1D(chpts, c(0.7, 0.2, 0.7), type = "geom"), normalize = FALSE),
    bern     = list(y = dataGenerator_1D(chpts, c(0.2, 0.8, 0.2), type = "bern"), normalize = FALSE),
    binom    = list(y = dataGenerator_1D(chpts, c(0.2, 0.8, 0.2), nbTrials = 10, type = "binom") / 10, normalize = FALSE),
    negbin   = list(y = dataGenerator_1D(chpts, c(0.7, 0.2, 0.7), nbSuccess = 10, type = "negbin") / 10, normalize = FALSE),
    variance = list(y = dataGenerator_1D(chpts, c(1, 5, 1), type = "variance"), normalize = FALSE)
  )

  for (model in names(cases)) {
    y <- cases[[model]]$y
    if (cases[[model]]$normalize) y <- data_normalization_1D(y, type = model)
    pen <- 2 * log(length(y))

    ref <- dust.1D(data = y, model = model, method = "PELT", penalty = pen, backend = "scalar")
    op  <- dust.1D(data = y, model = model, method = "OP", penalty = pen, backend = "scalar")

    expect_identical(op$changepoints, ref$changepoints, info = model)
    expect_equal(op$costQ, ref$costQ, tolerance = 1e-9, info = model)
    expect_true(all(op$nb == seq_along(op$nb)), info = model)
  }
})

test_that("dust.object.1D OP engine: one-shot, single-append and batched-append all agree", {
  set.seed(43)
  true_chpts <- c(300, 600, 900)
  y <- dataGenerator_1D(true_chpts, c(0, 2, -1), sdNoise = 1, type = "gauss")
  y <- data_normalization_1D(y, type = "gauss")
  pen <- 2 * log(length(y))

  # Batches deliberately don't align with true_chpts, to exercise resuming
  # the DP mid-segment, not just at a true change point.
  batches <- list(1:137, 138:700, 701:900)

  oneshot <- dust.1D(y, penalty = pen, model = "gauss", method = "OP", backend = "scalar")

  ob_single <- dust.object.1D(model = "gauss", method = "OP", backend = "scalar")
  ob_single$append_data(y, pen)
  ob_single$update_partition()
  res_single <- ob_single$get_partition()

  ob_batched <- dust.object.1D(model = "gauss", method = "OP", backend = "scalar")
  for (b in batches) {
    ob_batched$append_data(y[b], pen)
    ob_batched$update_partition()
  }
  res_batched <- ob_batched$get_partition()

  expect_identical(as.integer(res_single$changepoints), as.integer(oneshot$changepoints))
  expect_equal(res_single$costQ, oneshot$costQ, tolerance = 1e-9)
  expect_identical(as.integer(res_batched$changepoints), as.integer(oneshot$changepoints))
  expect_equal(res_batched$costQ, oneshot$costQ, tolerance = 1e-9)
  expect_true(all(unlist(res_batched$nb) == seq_along(unlist(res_batched$nb))))
})

test_that("an unrecognized method or model name errors instead of silently falling back", {
  expect_error(dust.object.1D(model = "gauss", method = "DUSTgs"))
  expect_error(dust.object.1D(model = "gauss", method = "det_DUST"))
  expect_error(dust.object.1D(model = "not_a_model"))
})

test_that("Variance model's statistic is the squared data, not the raw data", {
  # Regression guard for the one model whose `statistic()` differs from the
  # others (Gauss/Poisson/.../Negbin all record the raw datum).
  y <- c(1, -2, 3, -4)
  res <- dust.object.1D(model = "variance", backend = "scalar")
  res$append_data(y, 2 * log(length(y)))
  info <- res$get_info()
  expect_equal(diff(info$data_statistic), y^2)
})
