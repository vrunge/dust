# Independent O(n^2) reference, using direct sums rather than package policies.
ib_reference <- function(y, beta, model) {
  z <- if (model == "variance") y^2 else y
  f <- function(a) switch(model,
    gauss = a*a/2,
    poisson = if (a == 0) 0 else a*(log(a)-1),
    exp = -log(a)-1,
    geom = if (a == 1) 0 else (a-1)*log(a-1)-a*log(a),
    bern = if (a == 0 || a == 1) 0 else a*log(a)+(1-a)*log1p(-a),
    binom = if (a == 0 || a == 1) 0 else a*log(a)+(1-a)*log1p(-a),
    negbin = if (a == 0) 0 else a*log(a)-(1+a)*log1p(a),
    variance = -.5*(log(a)+1))
  q <- c(-beta, rep(NA_real_, length(y)))
  for (t in seq_along(y))
    q[t+1] <- min(vapply(0:(t-1), function(s)
      q[s+1]-(t-s)*f(mean(z[(s+1):t]))+beta, 0.))
  q[-1]
}

test_that("DUSTib scalar and Highway match independent unpruned optimization", {
  set.seed(1042)
  cases <- list(gauss=rnorm(45), poisson=rpois(45,.4), exp=rexp(45),
    geom=rgeom(45,.8)+1, bern=rbinom(45,1,.7), binom=rbinom(45,10,.7)/10,
    negbin=rnbinom(45,10,.8)/10, variance=rnorm(45,sd=2))
  for (model in names(cases)) for (beta in c(0,1,10)) {
    y <- cases[[model]]
    ref <- ib_reference(y,beta,model)
    for (backend in c("scalar", "highway")) {
      out <- dust.1D(y,beta,model,"DUSTib",backend=backend)
      expect_equal(as.numeric(out$costQ),ref,tolerance=2e-10,info=model)
    }
  }
})

test_that("DUSTib handles closed discrete boundaries and small positive statistics", {
  cases <- list(poisson=rep(c(0,0,1),11),geom=rep(c(1,1,2),11),
    bern=rep(c(0,1,1),11),binom=rep(c(0,.5,1),11),negbin=rep(c(0,0,.1),11),
    exp=rep(c(1e-20,2e-20,3e-20),11),variance=rep(c(1e-10,2e-10,3e-10),11))
  for (model in names(cases)) {
    y <- cases[[model]]; ref <- ib_reference(y,2,model)
    expect_equal(as.numeric(dust.1D(y,2,model,"DUSTib",backend="scalar")$costQ),ref,tolerance=2e-10)
    expect_equal(as.numeric(dust.1D(y,2,model,"DUSTib",backend = "highway")$costQ),ref,tolerance=2e-10)
    ob <- dust.object.1D(model,"DUSTib",backend = "highway")
    for (batch in list(1:1,2:14,15:length(y))) {
      ob$append_data(y[batch],2);ob$update_partition()
    }
    expect_equal(as.numeric(ob$get_partition()$costQ),ref,tolerance=2e-10)
  }
})

test_that("DUSTib rejects singular and invalid inputs before changing object state", {
  for (model in c("exp","variance")) for (backend in c("scalar", "highway")) {
    ob <- dust.object.1D(model,"DUSTib",backend=backend)
    ob$append_data(c(1,2),2);ob$update_partition()
    before <- ob$get_info()$data_length
    expect_error(ob$append_data(c(1,0),2),"strictly positive")
    expect_equal(ob$get_info()$data_length,before)
  }
  expect_error(dust.1D(c(0,2),1,"binom","DUSTib"),"proportions")
  expect_error(dust.1D(c(1,NA_real_),1,"gauss","DUSTib",backend = "highway"),"domain")
})

test_that("DUSTib Highway objects report the native backend when available", {
  ob <- dust.object.1D("gauss","DUSTib",backend = "highway")
  if (dust:::DUST.1D.HW.backend() == "highway") {
    expect_equal(ob$get_info()$backend,"highway")
    expect_equal(ob$get_info()$pruning_algo,"DUSTib")
  } else expect_s4_class(ob,"Rcpp_DUST_1D")
})
