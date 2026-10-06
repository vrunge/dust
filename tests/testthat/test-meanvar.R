reference_meanvar_cost <- function(y, penalty) {
  n <- length(y)
  sums <- c(0, cumsum(y))
  sums2 <- c(0, cumsum(y * y))
  q <- c(-penalty, rep(Inf, n))
  for (t in seq_len(n)) {
    for (s in 0:(t - 1L)) {
      len <- t - s
      if (len < 2L) next
      mean <- (sums[t + 1L] - sums[s + 1L]) / len
      variance <- (sums2[t + 1L] - sums2[s + 1L]) / len - mean * mean
      if (!(variance > 0)) next
      value <- q[s + 1L] + penalty + len * (1 + log(variance)) / 2
      if (value < q[t + 1L]) q[t + 1L] <- value
    }
  }
  q[-1L]
}

test_that("both meanVar methods and backends solve the same optimal partition", {
  set.seed(73)
  series <- list(
    c(rnorm(40), rnorm(45, 2, 1.5), rnorm(35, -1, 0.7)),
    rnorm(70),
    rep(c(-1, 1), 25)
  )
  for (y in series) {
    penalty <- 4 * log(length(y))
    expected <- reference_meanvar_cost(y, penalty)
    for (method in c("1D", "2D")) {
      scalar <- dust.meanVar(y, penalty, method, "scalar")
      highway <- dust.meanVar(y, penalty, method, "highway")
      expect_identical(names(scalar), c("changepoints", "lastIndexSet", "backend", "nb", "costQ"))
      expect_equal(scalar$costQ, expected, tolerance = 1e-6)
      expect_equal(highway$costQ, expected, tolerance = 1e-6)
      expect_identical(scalar$changepoints, highway$changepoints)
      expect_identical(scalar$lastIndexSet, highway$lastIndexSet)
      expect_identical(scalar$nb, highway$nb)
      expect_identical(scalar$backend, "scalar")
      expect_true(highway$backend %in% c("highway", "scalar"))
      expect_identical(highway$backend, dust:::DUST.1D.HW.backend())
    }
  }
})

test_that("meanVar objects resume from uneven batches for both backends", {
  set.seed(74)
  y <- c(rnorm(29), rnorm(40, 1, 2), rnorm(36, -1, 0.8))
  penalty <- 4 * log(length(y))
  for (method in c("1D", "2D")) for (backend in c("scalar", "highway")) {
    one_shot <- dust.meanVar(y, penalty, method, backend)
    object <- dust.object.meanVar(method, backend)
    object$append_data(y[1:29], penalty)
    object$update_partition()
    expect_length(object$get_partition()$costQ, 29)
    object$append_data(y[30:69], NULL)
    object$update_partition()
    object$append_data(y[70:105], NULL)
    object$update_partition()
    expect_equal(object$get_partition(), one_shot, tolerance = 1e-9)
    expect_identical(object$get_info()$backend, one_shot$backend)
  }
})

test_that("meanVar validates its public options and data", {
  expect_error(dust.meanVar(numeric()), "nonempty")
  expect_error(dust.meanVar(c(1, NA_real_)), "finite")
  expect_error(dust.meanVar(c(1, 2), penalty = -1), "penalty")
  expect_error(dust.meanVar(c(1, 2), method = "DUST1"), "arg")
  expect_error(dust.object.meanVar(backend = "auto"), "arg")
})

test_that("meanVar defaults to highway for both methods and object forms", {
  y <- c(-1, 1, -2, 2, -1, 1)
  for (method in c("1D", "2D")) {
    expected <- dust.meanVar(y, penalty = 2, method = method,
                             backend = "highway")
    expect_identical(dust.meanVar(y, penalty = 2, method = method), expected)
    expect_identical(dust.meanVar(y, penalty = 2, method = method,
                                  backend = "Highway"), expected)
    expect_identical(expected$backend, dust:::DUST.1D.HW.backend())
    object <- dust.object.meanVar(method = method)
    expect_identical(object$get_info()$backend, expected$backend)
    object$append_data(y, 2)
    object$update_partition()
    expect_identical(object$get_partition(), expected)
  }
})

test_that("meanVar reproduces the original dust pruning sets", {
  y <- sin((1:40) * 1.31) + cos((1:40) * 0.37)
  expected <- list(
    `1D` = list(
      nb = c(1,2,3,4,5,5,5,5,5,6,7,8,9,9,9,7,7,6,6,7,
             7,8,9,9,9,9,8,7,8,9,8,7,7,8,8,8,8,8,9,10),
      last = c(40,39,38,37,36,35,33,32,29,28,0)),
    `2D` = list(
      nb = c(1,2,3,4,5,5,5,5,5,6,6,7,8,8,8,7,6,6,6,7,
             7,8,8,7,6,7,7,7,8,9,7,6,7,8,6,7,7,8,9,10),
      last = c(40,39,38,37,36,35,33,32,29,28,0)))
  for (method in c("1D", "2D")) for (backend in c("scalar", "highway")) {
    result <- dust.meanVar(y, 4 * log(length(y)), method, backend)
    expect_identical(as.integer(result$nb), as.integer(expected[[method]]$nb))
    expect_identical(as.integer(result$lastIndexSet),
                     as.integer(expected[[method]]$last))
  }
})
