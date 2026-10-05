test_that("the two public dust interfaces select the requested backend", {
  y <- c(0, 0, 1, 1, 0, 0)
  penalty <- 1

  for (method in c("DUST", "DUSTib", "PELT", "OP")) {
    expect_identical(dust.1D(y, penalty, method = method),
                     dust.1D(y, penalty, method = method, backend = "highway"))
    expect_identical(dust.1D(y, penalty, method = method),
                     dust.1D(y, penalty, method = method, backend = "Highway"))
    scalar <- dust.object.1D(method = method, backend = "scalar")
    scalar$append_data(y, penalty)
    scalar$update_partition()
    expect_identical(scalar$get_info()$backend, "scalar")
    expect_identical(dust.1D(y, penalty, method = method, backend = "scalar"),
                     scalar$get_partition())

    ob <- dust.object.1D(method = method)
    ob$append_data(y, penalty)
    ob$update_partition()
    expect_identical(ob$get_info()$backend, dust:::DUST.1D.HW.backend())
    expect_equal(as.integer(ob$get_partition()$changepoints),
                 as.integer(dust.1D(y, penalty, method = method)$changepoints))
  }

  expect_error(dust.1D(y, penalty, backend = "auto"), "arg")
  expect_error(dust.object.1D(backend = "auto"), "arg")
})

test_that("dust.1D reports which backend actually ran, between lastIndexSet and nb", {
  y <- c(0, 0, 1, 1, 0, 0)
  penalty <- 1

  res_highway <- dust.1D(y, penalty, backend = "highway")
  expect_identical(names(res_highway)[3], "backend")
  expect_identical(res_highway$backend, dust:::DUST.1D.HW.backend())

  res_scalar <- dust.1D(y, penalty, backend = "scalar")
  expect_identical(names(res_scalar)[3], "backend")
  expect_identical(res_scalar$backend, "scalar")

  # Same field, same position, for every method -- including OP, which
  # (unlike DUST/PELT/DUSTib) has entirely separate scalar/Highway engines
  # (DUST_1D_OP_T / DUST_1D_HW_OP_T) with their own get_partition().
  for (method in c("DUST", "DUSTib", "PELT", "OP")) {
    for (backend in c("highway", "scalar")) {
      res <- dust.1D(y, penalty, method = method, backend = backend)
      expect_identical(names(res), c("changepoints", "lastIndexSet", "backend", "nb", "costQ"),
                        info = paste(method, backend))
    }
  }
})

test_that("only the two dust segmentation entry points are exported", {
  exports <- getNamespaceExports("dust")
  expect_true(all(c("dust.1D", "dust.object.1D") %in% exports))
  expect_false(any(c("DUST.1D.HW", "DUST.1D.HW.backend",
                     "dust.object.1D.HW") %in% exports))
})

test_that("Highway accepts all eight models with all four methods", {
  cases <- list(
    gauss = c(-1, 0, 1, 2, -1, 0),
    poisson = c(0, 1, 2, 0, 1, 2),
    exp = c(1, 2, 3, 1, 2, 3),
    geom = c(1, 2, 3, 1, 2, 3),
    bern = c(0, 1, 0, 1, 0, 1),
    binom = c(0, 0.5, 1, 0, 0.5, 1),
    negbin = c(0, 0.5, 1, 0, 0.5, 1),
    variance = c(1, 2, 3, 1, 2, 3)
  )
  for (model in names(cases)) for (method in c("DUST", "DUSTib", "PELT", "OP")) {
    y <- cases[[model]]
    hw <- dust.1D(y, penalty = 2, model = model, method = method,
                  backend = "highway")
    scalar <- dust.1D(y, penalty = 2, model = model, method = method,
                      backend = "scalar")
    expect_equal(hw$costQ, scalar$costQ, tolerance = 1e-9,
                 info = paste(model, method))
  }
})
