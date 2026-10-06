test_that("Highway and scalar backends match for every cost model and method", {
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

  # Both count models use their known size parameter to form normalized inputs.
  # Independent unpruned-DP correctness checks are in test-dustib.R.
  for (model in names(cases)) {
    y <- cases[[model]]$y
    if (cases[[model]]$normalize) y <- data_normalization_1D(y, type = model)
    pen <- 2 * log(length(y))
    for (method in c("DUST", "DUSTib", "PELT", "OP")) {
      ref <- dust.1D(data = y, model = model, method = method, penalty = pen, backend = "scalar")
      hw  <- dust.1D(y, pen, model = model, method = method, backend = "highway")
      info <- paste(model, method)
      expect_identical(as.integer(hw$changepoints), as.integer(ref$changepoints), info = info)
      expect_equal(hw$costQ, ref$costQ, tolerance = 1e-9, info = info)
    }
  }
})

test_that("DUST.1D.HW OP engine never prunes and matches PELT's costQ exactly, for every cost model", {
  # Regression guard for DUST_1D_HW_OP_T / run_OP_HW (hw_op_step in
  # DUST_1D_HW.cpp): the old HW method="OP" (methodCode==2) still ran
  # smallest_index_prune+compact every step, so this checks both
  # ground-truth agreement with PELT and that nb is trivially t at every
  # step (the active set is never reduced).
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
    op  <- dust.1D(y, pen, model = model, method = "OP", backend = "highway")

    expect_identical(as.integer(op$changepoints), as.integer(ref$changepoints), info = model)
    expect_equal(op$costQ, ref$costQ, tolerance = 1e-9, info = model)
    expect_true(all(unlist(op$nb) == seq_along(unlist(op$nb))), info = model)
  }
})

test_that("dust.object.1D.HW: one-shot, single-append and uneven batched-append all agree", {
  set.seed(43)
  true_chpts <- c(300, 600, 900)
  y <- dataGenerator_1D(true_chpts, c(0, 2, -1), sdNoise = 1, type = "gauss")
  y <- data_normalization_1D(y, type = "gauss")
  pen <- 2 * log(length(y))

  # Batches deliberately don't align with true_chpts, to exercise resuming
  # the DP mid-segment, not just at a true change point.
  batches <- list(1:137, 138:700, 701:900)

  for (method in c("DUST", "DUSTib", "PELT", "OP")) {
    oneshot <- dust.1D(y, pen, model = "gauss", method = method, backend = "highway")

    ob_single <- dust.object.1D(model = "gauss", method = method, backend = "highway")
    ob_single$append_data(y, pen)
    ob_single$update_partition()
    res_single <- ob_single$get_partition()

    ob_batched <- dust.object.1D(model = "gauss", method = method, backend = "highway")
    for (b in batches) {
      ob_batched$append_data(y[b], pen)
      ob_batched$update_partition()
    }
    res_batched <- ob_batched$get_partition()

    info <- method
    expect_identical(as.integer(res_single$changepoints), as.integer(oneshot$changepoints), info = info)
    expect_equal(res_single$costQ, oneshot$costQ, tolerance = 1e-9, info = info)
    expect_identical(as.integer(res_batched$changepoints), as.integer(oneshot$changepoints), info = info)
    expect_equal(res_batched$costQ, oneshot$costQ, tolerance = 1e-9, info = info)
  }
})

test_that("DUST.1D.HW.backend reports a valid backend", {
  expect_true(dust:::DUST.1D.HW.backend() %in% c("highway", "scalar"))
})

test_that("dust.object.1D.HW uses the available backend for DUSTib", {
  ob <- dust.object.1D(model = "gauss", method = "DUSTib", backend = "highway")
  expect_s4_class(ob, if (dust:::DUST.1D.HW.backend() == "highway") "Rcpp_DUST_1D_HW_Obj" else "Rcpp_DUST_1D")
})

test_that("DUST.1D.HW errors on an unrecognized model or method instead of silently falling back", {
  y <- data_normalization_1D(dataGenerator_1D(c(50, 100), c(0, 1), type = "gauss"), type = "gauss")
  expect_error(dust.1D(y, model = "not_a_model", backend = "highway"))
  expect_error(dust.1D(y, method = "not_a_method", backend = "highway"))
})
