# Optional solver-level checks against analytic maxima, with no public test API.
# Run from the package root: Rscript tools/check-md-search.R [package-directory]
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
source_file <- normalizePath(file.path(root, "src", "MD_DUST.cpp"),
                             winslash = "/", mustWork = TRUE)
code <- paste0('
// [[Rcpp::plugins(cpp17)]]
#include ', encodeString(source_file, quote = '"'), '
template<class Model>
Rcpp::List probe_model(Rcpp::NumericVector S, Rcpp::NumericMatrix M,
                       double Q, Rcpp::NumericVector U,
                       Rcpp::NumericVector x, int loops, bool qn) {
  Decision<Model> d;
  d.dimension = S.size(); d.constraints = U.size(); d.highway = false;
  d.a = Rcpp::as<std::vector<double>>(S);
  d.matrix.assign(M.begin(), M.end());
  d.u = Rcpp::as<std::vector<double>>(U); d.c = Q;
  d.scratch.resize(d.dimension);
  double scale = 0;
  const double value = d.value(Rcpp::as<std::vector<double>>(x), scale);
  return Rcpp::List::create(Rcpp::Named("prune") = iterative_search(d, loops, qn),
                           Rcpp::Named("value") = value);
}
// [[Rcpp::export]]
Rcpp::List probe_md(int model, Rcpp::NumericVector S, Rcpp::NumericMatrix M,
                   double Q, Rcpp::NumericVector U, Rcpp::NumericVector x,
                   int loops, bool qn) {
  if (M.nrow() != S.size() || M.ncol() != U.size() || x.size() != U.size())
    Rcpp::stop("incompatible dimensions");
  switch(model) {
    case 0: return probe_model<GaussPolicy>(S,M,Q,U,x,loops,qn);
    case 1: return probe_model<PoissonPolicy>(S,M,Q,U,x,loops,qn);
    case 2: return probe_model<ExpPolicy>(S,M,Q,U,x,loops,qn);
    case 3: return probe_model<GeomPolicy>(S,M,Q,U,x,loops,qn);
    case 4: return probe_model<BernPolicy>(S,M,Q,U,x,loops,qn);
    case 5: return probe_model<BinomPolicy>(S,M,Q,U,x,loops,qn);
    case 6: return probe_model<NegbinPolicy>(S,M,Q,U,x,loops,qn);
    case 7: return probe_model<VariancePolicy>(S,M,Q,U,x,loops,qn);
  }
  Rcpp::stop("invalid model");
}
')
Rcpp::sourceCpp(code = code, showOutput = FALSE)

models <- c("gauss", "poisson", "exp", "geom", "bern", "binom",
            "negbin", "variance")
conjugate <- function(model, z) switch(model,
  gauss = z^2 / 2, poisson = z * (log(z) - 1), exp = -log(z) - 1,
  geom = (z - 1) * log1p(-1 / z) - log(z),
  bern =, binom = z * log(z) + (1 - z) * log1p(-z),
  negbin = z * (log(z) - log1p(z)) - log1p(z),
  variance = -(log(z) + 1) / 2)
theta <- function(model, z) switch(model,
  gauss = z, poisson = log(z), exp = -1 / z, geom = log1p(-1 / z),
  bern =, binom = log(z) - log1p(-z),
  negbin = log(z) - log1p(z), variance = -0.5 / z)

checks <- 0L
check <- function(model, S, M, Q, U, x, positive, value = NULL) {
  for (qn in c(FALSE, TRUE)) {
    result <- probe_md(match(model, models) - 1L, S, M, Q, U, x, 400L, qn)
    stopifnot(identical(result$prune, positive))
    if (!is.null(value)) stopifnot(isTRUE(all.equal(result$value, value,
                                                  tolerance = 1e-11)))
    checks <<- checks + 1L
  }
}

# Stationarity at x*=1 and concavity give an analytic global maximum.
# Include decreasing means (which can send trial points outside the domain),
# independent constraints, and repeated columns with singular curvature.
for (model in models) for (increasing in c(FALSE, TRUE)) {
  bounded <- model %in% c("bern", "binom")
  S <- if (bounded) 0.6 else if (model == "geom") 3 else 2
  z <- if (bounded) 0.25 else if (model == "geom") 1.5 else 0.7
  if (increasing) {
    S <- if (bounded) 0.4 else 3
    z <- if (bounded) 0.65 else 4
  }
  m <- z - S
  u <- -m * theta(model, z)
  peak <- -conjugate(model, z) - u
  gain <- peak + conjugate(model, S)
  stopifnot(gain > 0)
  for (margin in c(-gain / 4, 0, gain / 4)) {
    check(model, S, matrix(m), peak - margin, u, 1, margin > 0, margin)
    check(model, rep(S, 2), diag(rep(m, 2)), 2 * peak - margin,
          rep(u, 2), c(1, 1), margin > 0, margin)
    check(model, S, matrix(c(m, m), 1), peak - margin,
          rep(u, 2), c(0.5, 0.5), margin > 0, margin)
  }
  # Even a negative multiplier with a positive extrapolated value is invalid.
  check(model, S, matrix(m), peak - gain / 4, u, -1, TRUE, -Inf)
  if (!increasing && model != "gauss")
    check(model, S, matrix(m), peak - gain / 4, u, 2, TRUE, -Inf)

  if (!increasing && model %in% c("poisson", "geom", "bern", "binom", "negbin")) {
    endpoint <- if (model == "geom") 1 else 0
    # x1 must be zero; the second coordinate must still find the witness.
    M <- matrix(c(-1, 0, 0, m), 2)
    check(model, c(endpoint, S), M, peak - gain / 4,
          c(-100, u), c(1, 0), TRUE, -Inf)
    check(model, endpoint, matrix(-1), 1, -100, 0, FALSE, -1)
    if (bounded) {
      M[1, 1] <- 1
      check(model, c(1, S), M, peak - gain / 4,
            c(-100, u), c(1, 0), TRUE, -Inf)
    }
  }
}

# Independent Gaussian oracle: enumerate all faces of the orthant.
set.seed(196)
for (p in 2:4) for (repetition in 1:8) {
  M <- diag(p) + matrix(rnorm(p * p, sd = 0.2), p)
  S <- rnorm(p)
  H <- crossprod(M)
  h <- rnorm(p)
  candidates <- list(rep(0, p))
  for (mask in seq_len(2^p - 1L)) {
    I <- which(as.logical(intToBits(mask)[seq_len(p)]))
    x <- numeric(p)
    x[I] <- solve(H[I, I, drop = FALSE], h[I])
    if (all(x >= 0)) candidates[[length(candidates) + 1L]] <- x
  }
  gains <- vapply(candidates, function(x) sum(h * x) - sum(x * (H %*% x)) / 2, 0.0)
  best <- which.max(gains)
  gain <- gains[best]
  U <- as.vector(-crossprod(M, S) - h)
  for (margin in c(-0.01, 0, 0.2 * gain + 0.01))
    check("gauss", S, M, -sum(S^2) / 2 + gain - margin, U,
          candidates[[best]], margin > 0, margin)
}

# Flat conjugate contribution: an unbounded positive direction, and a
# stationary negative decision that must be retained.
check("gauss", 0, matrix(0), 1, -1, 0, TRUE, -1)
check("gauss", 0, matrix(0), 1, 0, 0, FALSE, -1)
cat(checks, "analytic decision-search checks passed.\n")
