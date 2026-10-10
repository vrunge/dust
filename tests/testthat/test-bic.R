test_that("the default penalties are the BIC", {
  set.seed(1)
  n <- 200
  one <- dust.object.1D()
  one$append_data(rnorm(n), NULL)
  expect_equal(one$get_info()$current_penalty, 2 * log(n))
  md <- dust.object.MD()
  md$append_data(matrix(rnorm(3 * n), nrow = 3), NULL)
  expect_equal(md$get_info()$current_penalty, 4 * log(n))
  mv <- dust.object.meanVar()
  mv$append_data(rnorm(n), NULL)
  expect_equal(mv$get_info()$current_penalty, 3 * log(n))
})

test_that("binom and negbin: size gives the penalty and the costs on the -2 log-likelihood scale", {
  set.seed(2)
  n <- 240
  size <- 7
  y <- dataGenerator_1D(chpts = c(80, 160, n), parameters = c(0.6, 0.2, 0.5), nbTrials = size,
                        nbSuccess = size, type = "binom")
  z <- dataGenerator_1D(chpts = c(80, 160, n), parameters = c(0.6, 0.2, 0.5), nbSuccess = size, type = "negbin")
  for (case in list(list("binom", y), list("negbin", z)))
  {
    model <- case[[1]]
    counts <- case[[2]]
    res <- dust.1D(counts, model = model, size = size)
    manual <- dust.1D(counts / size, 2 * log(n) / size, model = model)
    expect_equal(res$changepoints, manual$changepoints, info = model)
    expect_equal(res$costQ, manual$costQ * size, info = model)
    expect_equal(dust.object.1D(model, size = size)$dust(counts, NULL)$costQ, res$costQ, info = model)

    # costQ = -2 log-likelihood (up to the binomial coefficients) + penalties
    ends <- c(0, res$changepoints)
    loglik <- 0
    for (i in seq_len(length(ends) - 1L))
    {
      seg <- counts[(ends[i] + 1L):ends[i + 1L]]
      m <- mean(seg)
      loglik <- loglik + if (model == "binom")
        sum(dbinom(seg, size, m / size, log = TRUE)) - sum(lchoose(size, seg))
      else sum(dnbinom(seg, size, mu = m, log = TRUE)) -
        sum(lgamma(seg + size) - lgamma(size) - lgamma(seg + 1))
    }
    k <- length(res$changepoints) - 1L
    expect_equal(tail(res$costQ, 1), -2 * loglik + k * 2 * log(n), info = model)
  }
})

test_that("size is checked", {
  expect_error(dust.1D(rpois(50, 3), model = "poisson", size = 10))
  expect_error(dust.1D(c(1, 2, 3), model = "binom", size = 2))
  expect_error(dust.1D(c(1, 2, 3), model = "binom", size = -1))
  expect_error(dust.1D(c(0.5, 1, 2), model = "negbin", size = 3))
})
