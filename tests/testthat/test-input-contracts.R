test_that("an empty first append does not initialize or corrupt a 1D object", {
  y <- c(0, 1, 0, 2, 1)
  for (backend in c("scalar", "highway"))
    for (method in c("DUST", "DUSTib", "PELT", "OP")) {
      object <- dust.object.1D(method = method, backend = backend)
      object$update_partition()
      expect_error(object$get_partition(), "append data")
      object$append_data(numeric(), NULL)
      expect_error(object$get_partition(), "append data")
      object$append_data(y, 2)
      expect_error(object$get_partition(), "update_partition")
      object$update_partition()
      expect_identical(object$get_partition(),
                       dust.1D(y, 2, method = method, backend = backend),
                       info = paste(backend, method))
      expect_identical(names(object$get_info()),
                       c("backend", "data_statistic", "data_length",
                         "current_penalty", "model", "pruning_algo"))
    }
})

test_that("all 1D methods enforce the observation domain and penalty", {
  invalid <- list(poisson = c(1, -1), exp = c(1, 0),
                  geom = c(1, 0), bern = c(0, 2), binom = c(0, 2),
                  negbin = c(0, -1), variance = c(1, 0),
                  gauss = c(1, NA_real_))
  for (backend in c("scalar", "highway"))
    for (method in c("DUST", "DUSTib", "PELT", "OP")) {
      for (model in names(invalid))
        expect_error(dust.1D(invalid[[model]], 2, model, method, backend),
                     "domain|finite", info = paste(backend, method, model))
      for (penalty in c(-1, Inf, NA_real_))
        expect_error(dust.1D(c(0, 1), penalty, backend = backend,
                             method = method), "penalty")
      object <- dust.object.1D(method = method, backend = backend)
      object$append_data(c(0, 1), 2)
      expect_error(object$append_data(c(1, 2), 3), "cannot change")
      expect_equal(object$get_info()$data_length, 2)
    }
})

test_that("meanVar online penalty is fixed and degenerate output is explicit", {
  for (backend in c("scalar", "highway")) {
    object <- dust.object.meanVar(backend = backend)
    object$append_data(c(0, 1), 2)
    expect_error(object$get_partition(), "update_partition")
    object$update_partition()
    expect_error(object$append_data(c(2, 3), 3), "cannot change")
    expect_equal(object$get_info()$data_length, 2)
  }
  expect_error(dust.meanVar(1, 2), "no finite")
  expect_error(dust.meanVar(rep(1, 4), 2), "no finite")
})
