# Solver-level checks for the Gaussian joint maximum. Run from the package
# source directory with: Rscript tools/check-md-exact.R
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", mustWork = TRUE)
source_file <- normalizePath(file.path(root, "src", "MD_DUST.cpp"),
                             winslash = "/", mustWork = TRUE)
code <- paste0('
// [[Rcpp::plugins(cpp17)]]
#include ', encodeString(source_file, quote = '"'), '
using namespace dust_md;
// [[Rcpp::export]]
Rcpp::List probe_gauss_joint(Rcpp::NumericVector S, Rcpp::NumericMatrix M,
                             double Q, Rcpp::NumericVector U) {
  if (M.nrow() != S.size() || M.ncol() != U.size())
    Rcpp::stop("incompatible dimensions");
  Decision<GaussPolicy> d;
  d.dimension = S.size(); d.constraints = U.size(); d.highway = false;
  d.a = Rcpp::as<std::vector<double>>(S);
  d.matrix.assign(M.begin(), M.end());
  d.u = Rcpp::as<std::vector<double>>(U); d.c = Q;
  d.scratch.resize(d.dimension);
  const auto result = gauss_joint_maximum(d);
  return Rcpp::List::create(
    Rcpp::Named("kind") = static_cast<int>(result.kind),
    Rcpp::Named("point") = result.point,
    Rcpp::Named("value") = result.value,
    Rcpp::Named("prune") = gauss_exact_search(d));
}
')
Rcpp::sourceCpp(code = code, showOutput = FALSE)

# A separate R implementation enumerates the KKT faces for regular cases.
oracle <- function(S, M, Q, U) {
  p <- ncol(M)
  G <- crossprod(M)
  h <- as.vector(-crossprod(M, S) - U)
  solutions <- list()
  for (mask in 0:(2^p - 1L)) {
    I <- which(as.logical(intToBits(mask)[seq_len(p)]))
    x <- numeric(p)
    if (length(I)) x[I] <- solve(G[I, I, drop = FALSE], h[I])
    if (any(x < 0)) next
    g <- as.vector(h - G %*% x)
    if (length(I) && any(abs(g[I]) > 1e-9)) next
    inactive <- setdiff(seq_len(p), I)
    if (length(inactive) && any(g[inactive] > 1e-9)) next
    solutions[[length(solutions) + 1L]] <- x
  }
  stopifnot(length(solutions) > 0L)
  x <- solutions[[1L]]
  list(point = x, value = -sum((S + M %*% x)^2) / 2 - Q - sum(U * x))
}

set.seed(2791)
checks <- 0L
for (p in 1:10) for (rep in seq_len(if (p <= 5) 60L else 8L)) {
  M <- diag(p) + matrix(rnorm(p * p, sd = 0.3), p)
  S <- rnorm(p)
  U <- rnorm(p)
  Q <- rnorm(1)
  expected <- oracle(S, M, Q, U)
  actual <- probe_gauss_joint(S, M, Q, U)
  if (actual$kind != 0L || !all(actual$point >= 0) ||
      !isTRUE(all.equal(actual$value, expected$value, tolerance = 1e-8)))
    stop(sprintf("p=%d replicate=%d: S=%s M=%s Q=%s U=%s expected %.16g at %s, got kind %d, value %.16g at %s",
                 p, rep, paste(S, collapse = ","), paste(as.vector(M), collapse = ","),
                 Q, paste(U, collapse = ","),
                 expected$value, paste(expected$point, collapse = ","),
                 actual$kind, actual$value, paste(actual$point, collapse = ",")))
  # Values near zero are checked separately, since the package uses a
  # scale-aware floating-point positivity guard.
  if (abs(expected$value) > 1e-6)
    stopifnot(identical(actual$prune, expected$value > 0))
  checks <- checks + 1L
}

# The joint maximum is positive even though both axis maxima are negative.
M <- cbind(c(1, 0), c(-0.5, sqrt(0.75)))
actual <- probe_gauss_joint(c(0, 0), M, 1, c(-1, -1))
stopifnot(actual$kind == 0L, actual$prune,
          isTRUE(all.equal(actual$value, 1, tolerance = 1e-12)),
          all(abs(actual$point - c(2, 2)) < 1e-10))
checks <- checks + 1L

# The unconstrained critical point has a negative coordinate. The true
# constrained solution is on a different face than componentwise clipping.
M <- cbind(c(sqrt(2), 0), c(1 / sqrt(2), sqrt(1.5)))
actual <- probe_gauss_joint(c(0, 0), M, 0, c(0, -1))
stopifnot(actual$kind == 0L,
          all(abs(actual$point - c(0, 0.5)) < 1e-10),
          isTRUE(all.equal(actual$value, 0.25, tolerance = 1e-12)))
checks <- checks + 1L

# Dependent columns: a finite face maximum, a mixed-sign null direction
# with zero slope, and two genuinely unbounded decisions.
cases <- list(
  list(S = 0, M = matrix(c(1, 1), 1), Q = 0.1, U = c(-1, -1),
       kind = 0L, value = 0.4, prune = TRUE),
  list(S = 0, M = matrix(c(1, -1), 1), Q = 0.3, U = c(-1, 1),
       kind = 0L, value = 0.2, prune = TRUE),
  list(S = 0, M = matrix(c(1, -1), 1), Q = 1, U = c(-1, -1),
       kind = 1L, value = Inf, prune = TRUE),
  list(S = 0, M = matrix(0, 1), Q = 1, U = -1,
       kind = 1L, value = Inf, prune = TRUE),
  list(S = c(0, 0), M = cbind(c(1, 0), c(0, 1), c(-1, -1)),
       Q = 1, U = rep(-1, 3), kind = 1L, value = Inf, prune = TRUE),
  list(S = 0, M = matrix(0, 1), Q = 1, U = 1,
       kind = 0L, value = -1, prune = FALSE))
for (case in cases) {
  actual <- probe_gauss_joint(case$S, case$M, case$Q, case$U)
  stopifnot(actual$kind == case$kind,
            identical(actual$prune, case$prune),
            all(actual$point >= 0))
  if (is.finite(case$value))
    stopifnot(isTRUE(all.equal(actual$value, case$value, tolerance = 1e-10)))
  else stopifnot(actual$value > 0)
  checks <- checks + 1L
}
cat(checks, "Gaussian joint-maximum checks passed.\n")
