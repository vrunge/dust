
###
### Section 1. INTRODUCTION (Poisson model)
###

#########
# system time, limit = 10s
# Three experiments: 10, 100 and 1000 equal-length segments (9, 99, 999
# true changepoints), rates alternating 2/8. For each, find the largest
# n every algorithm can segment within a 10s budget.
#########
# BS    (changepoint::cpt.meanvar, method="BinSeg")
# OP    best of: changepoint::cpt.meanvar(method="SegNeigh") vs
#       dust::dust.1D(model="poisson", method="OP", backend = "highway")
# PELT  best of: changepoint::cpt.meanvar(method="PELT") vs
#       dust::dust.1D(model="poisson", method="PELT", backend = "highway")
# GFPOP (gfpop::gfpop -- the Poisson analogue of FPOP)
# DUST  (dust::dust.1D, model="poisson", method="DUST", backend = "highway")
#
# OP and PELT are no longer single-implementation rows: dust's own OP
# and PELT (method="OP"/"PELT" on dust.1D with backend = "highway", generalized to all 8
# models and validated against the scalar backend across a 128-combination
# stress matrix earlier this session) are real competitors in their own
# right, not just dust's headline DUST algorithm. Both implementations
# are run and bisected independently; the larger max_n wins and is the
# one reported, with which implementation won recorded alongside it.

library(gfpop)
library(changepoint)
library(dust)
library(callr)

# Re-run on Apple M4 Pro (14-core: 10 performance + 4 efficiency), 24GB,
# macOS 26.6.2 (see Section1_introduction_gauss.R for the same note).
#
# Unlike the Gaussian experiment, OP and PELT here include both dust
# and the external reference implementations (changepoint::cpt.meanvar).
# gfpop (graph-
# constrained FPOP) stands in for plain FPOP, which is Gaussian-only.
#
# SegNeigh is the classical Segment Neighbourhood algorithm: O(Q * n^2)
# in changepoint's implementation (full DP search over every number of
# segments up to Q). This means its 10s-budget capacity can be orders of
# magnitude smaller than every other algorithm here, especially once Q
# (= nb_seg) grows -- the bisection brackets below start deliberately
# small for SegNeigh specifically, to avoid a single wildly-overshot
# evaluation costing minutes (exactly what happened once in the
# Gaussian run's OP bracket-finding step with a too-generous starting
# n_hi: worth remembering as a general lesson, not just a SegNeigh one).
#
# DUST uses dust.1D(..., backend = "highway") on data normalized the way dust expects
# (data_normalization_1D(y, type="poisson")); every competitor here
# gets the raw (unnormalized) Poisson counts, since none of them share
# dust's normalization convention. Same BIC-type penalty, 2*log(n),
# for all five either way.

################################################################################
##### Generic bisection search for "max n within budget seconds"
################################################################################

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  rates <- rep(c(2, 8), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = rates, type = "poisson")
}

set.seed(2026)

# changepoint's recursive BinSeg and SegNeigh (Poisson) can hit a C stack
# overflow for large n: confirmed by direct reproduction (n=5e6, nb_seg=10
# crashes reliably regardless of seed), not just a rare/borderline case
# as first suspected. Worse, the crash bypasses ordinary tryCatch: R's
# "C stack usage too close to the limit" protection can trigger too late
# for its own error-raising machinery to safely unwind through registered
# handlers, so a plain tryCatch in this same process cannot be trusted to
# catch it (confirmed: it didn't, in-process, even with a hard-limit-sized
# OS stack via ulimit -s).
#
# The robust fix is process isolation, not a bigger stack: run BinSeg and
# SegNeigh in a disposable callr subprocess (see run_BS/run_OP below). If
# that subprocess's C stack overflows, only it dies -- callr::r() detects
# this from the healthy parent and returns NULL (confirmed empirically:
# not a catchable error in the parent, just a NULL return), which is
# trivial to treat as "this n is not usable" (Inf elapsed, i.e. certainly
# over budget) without risking the whole multi-hour search.
time_it <- function(run_fun, n, nb_seg) {
  data <- build_data(n, nb_seg)
  t <- tryCatch(run_fun(data, n, nb_seg), error = function(e) Inf)
  if (is.null(t)) Inf else t
}

# n_floor: build_data() needs n >= nb_seg (one point per segment at an
# absolute minimum) to produce valid, strictly-positive changepoints --
# the shrink loop below must never cross that floor, however slow/Inf
# the measurements near it are.
find_max_n <- function(run_fun, nb_seg, n_lo, n_hi, label, budget = 10, tol = 0.02) {
  n_floor <- 2 * nb_seg
  t_lo <- time_it(run_fun, n_lo, nb_seg)
  while (t_lo > budget) {
    if (n_lo <= n_floor) { cat(sprintf("[%s | nb_seg=%d] hit n_floor=%.0f still over budget (%.2fs) -- reporting floor\n", label, nb_seg, n_floor, t_lo)); return(n_floor) }
    n_hi <- n_lo; n_lo <- max(n_floor, n_lo / 2); t_lo <- time_it(run_fun, n_lo, nb_seg)
  }
  t_hi <- time_it(run_fun, n_hi, nb_seg)
  while (t_hi < budget) { n_lo <- n_hi; n_hi <- n_hi * 1.5; t_hi <- time_it(run_fun, n_hi, nb_seg) }
  cat(sprintf("[%s | nb_seg=%d] bracket: n=%.0f (%.2fs) .. n=%.0f (%.2fs)\n", label, nb_seg, n_lo, t_lo, n_hi, t_hi))
  while ((n_hi - n_lo) / n_hi > tol) {
    n_mid <- (n_lo + n_hi) / 2
    t_mid <- time_it(run_fun, n_mid, nb_seg)
    cat(sprintf("[%s | nb_seg=%d] n=%.0f -> %.3fs\n", label, nb_seg, n_mid, t_mid))
    if (t_mid <= budget) n_lo <- n_mid else n_hi <- n_mid
  }
  cat(sprintf("[%s | nb_seg=%d] FINAL max n within %ds: %.0f\n\n", label, nb_seg, budget, n_lo))
  n_lo
}

# `baseline` runs first, from its own (hand-calibrated) bracket. Every
# other candidate is then warm-started AT the baseline's optimum --
# find_max_n's own grow-loop extends upward from there "if necessary"
# (i.e. if this candidate beats the baseline), and its shrink-loop
# handles the opposite case (this candidate slower than the baseline at
# that n) exactly as it would for any ordinary over-estimate. This
# rescales automatically with nb_seg (no stale fixed brackets), answers
# the actual question -- "is this candidate better than the baseline?"
# -- directly, and costs nothing in safety: unlike changepoint's
# BinSeg/SegNeigh, dust has no crash risk at any n, so anchoring right
# at the baseline's number is fine in either direction.
best_of <- function(candidates, nb_seg, label, baseline) {
  order_names <- c(baseline, setdiff(names(candidates), baseline))
  best_n <- -Inf; best_name <- NA_character_; baseline_n <- NULL
  for (nm in order_names) {
    cnd <- candidates[[nm]]
    n_lo <- cnd$n_lo; n_hi <- cnd$n_hi
    if (!is.null(baseline_n)) {
      n_lo <- baseline_n
      n_hi <- max(cnd$n_hi, baseline_n * 20)
    }
    n <- find_max_n(cnd$run_fun, nb_seg, n_lo = n_lo, n_hi = n_hi, label = paste0(label, "/", nm))
    cat(sprintf("  [%s/%s] max_n = %.0f\n", label, nm, n))
    if (nm == baseline) baseline_n <- n
    if (n > best_n) { best_n <- n; best_name <- nm }
  }
  cat(sprintf("[%s] BEST: %s (max_n=%.0f)\n\n", label, best_name, best_n))
  list(max_n = best_n, impl = best_name)
}

# Every run_* returns ELAPSED SECONDS directly (not a result object) --
# BS/OP time themselves inside their own subprocess (excluding that
# subprocess's own R startup, which would otherwise inflate the
# measurement by its ~0.3-0.5s launch overhead on every single call) and
# report NULL back out if they crashed; time_it() (above) turns that NULL
# into Inf. Q well above nb_seg (BIC can find a handful more changepoints
# than the true count) but capped well below n: Q approaching or
# exceeding n -- more candidate changepoints than data points -- is a
# second, independent crash trigger from the large-n one, hit directly
# at tiny n (Q=30 against n=20-30 crashed SegNeigh immediately).
run_BS <- function(data, n, nb_seg) {
  callr::r(function(data, n, nb_seg) {
    library(changepoint)
    Q <- min(2*nb_seg + 10, max(nb_seg, floor(n/3)))
    t0 <- Sys.time()
    changepoint::cpt.meanvar(data, test.stat = "Poisson", method = "BinSeg", penalty = "BIC", Q = Q)
    as.numeric(Sys.time() - t0, units = "secs")
  }, args = list(data = data, n = n, nb_seg = nb_seg))
}
run_OP_cp <- function(data, n, nb_seg) {
  callr::r(function(data, n, nb_seg) {
    library(changepoint)
    Q <- min(2*nb_seg + 10, max(nb_seg, floor(n/3)))
    t0 <- Sys.time()
    changepoint::cpt.meanvar(data, test.stat = "Poisson", method = "SegNeigh", penalty = "BIC", Q = Q)
    as.numeric(Sys.time() - t0, units = "secs")
  }, args = list(data = data, n = n, nb_seg = nb_seg))
}
run_PELT_cp <- function(data, n, nb_seg) {
  system.time(cpt.meanvar(data, test.stat = "Poisson", method = "PELT", penalty = "BIC"))[["elapsed"]]
}
# dust's own OP/PELT (dust.1D with backend = "highway", generalized to all 8 models and all 4
# methods earlier this session, validated against the scalar backend across a
# 128-combination stress matrix): real competitors here, not just the
# headline DUST row. Same normalization convention as run_DUST.
run_OP_dust <- function(data, n, nb_seg) {
  yn <- data_normalization_1D(data, type = "poisson")
  system.time(dust::dust.1D(yn, 2 * log(n), model = "poisson", method = "OP", backend = "highway"))[["elapsed"]]
}
run_PELT_dust <- function(data, n, nb_seg) {
  yn <- data_normalization_1D(data, type = "poisson")
  system.time(dust::dust.1D(yn, 2 * log(n), model = "poisson", method = "PELT", backend = "highway"))[["elapsed"]]
}
run_GFPOP <- function(data, n, nb_seg) {
  g <- graph(type = "std", penalty = 2 * log(n))
  system.time(gfpop(data, mygraph = g, type = "poisson"))[["elapsed"]]
}
run_DUST <- function(data, n, nb_seg) {
  yn <- data_normalization_1D(data, type = "poisson")
  system.time(dust::dust.1D(yn, 2 * log(n), model = "poisson", method = "DUST", backend = "highway"))[["elapsed"]]
}

results <- data.frame(nb_seg = integer(0), algo = character(0), max_n = numeric(0), impl = character(0))

for (nb_seg in c(10, 100, 1000)) {
  n_BS <- find_max_n(run_BS, nb_seg, n_lo = 1e4, n_hi = 5e6, label = "BS")

  # dust's OP/PELT are NOT DUST-scale here: OP by construction does
  # almost no pruning (only the unconditional smallest-index check), and
  # for this Poisson rate pattern the active set grows roughly linearly
  # with n (measured: final nb ~ 0.1*n at n=5e4..2e5), i.e. close to
  # O(n^2) overall -- capacity is of the same order as the changepoint
  # versions, not anywhere near DUST's 1e7-1e8 scale. changepoint runs
  # first (baseline) from its own hand-calibrated bracket; dust is then
  # warm-started at the baseline's optimum and grows/shrinks from there
  # -- see best_of()'s comment for why this beats giving dust its own
  # independent guessed bracket.
  op_best <- best_of(list(
    changepoint = list(run_fun = run_OP_cp,    n_lo = max(2e1, 3*nb_seg), n_hi = max(2e3, 300*nb_seg)),
    dust       = list(run_fun = run_OP_dust, n_lo = 1e4,                n_hi = 1e6)
  ), nb_seg, "OP", baseline = "changepoint")

  pelt_best <- best_of(list(
    changepoint = list(run_fun = run_PELT_cp,    n_lo = 1e4, n_hi = 5e6),
    dust       = list(run_fun = run_PELT_dust, n_lo = 1e4, n_hi = 5e5)
  ), nb_seg, "PELT", baseline = "changepoint")

  n_GFPOP <- find_max_n(run_GFPOP, nb_seg, n_lo = 1e5,  n_hi = 2e7,  label = "GFPOP")
  n_DUST  <- find_max_n(run_DUST,  nb_seg, n_lo = 1e6,  n_hi = 1e8,  label = "DUST")

  results <- rbind(results, data.frame(
    nb_seg = nb_seg,
    algo   = c("BS", "OP", "PELT", "GFPOP", "DUST"),
    max_n  = c(n_BS, op_best$max_n, pelt_best$max_n, n_GFPOP, n_DUST),
    impl   = c("changepoint", op_best$impl, pelt_best$impl, "gfpop", "dust")
  ))
}

################################################################################
cat("================ SUMMARY (max n within 10s) ================\n")
print(results, row.names = FALSE)

###
### nb_seg =   10 : BS = ___ , OP = ___ , PELT = ___ , GFPOP = ___ , DUST = ___
### nb_seg =  100 : BS = ___ , OP = ___ , PELT = ___ , GFPOP = ___ , DUST = ___
### nb_seg = 1000 : BS = ___ , OP = ___ , PELT = ___ , GFPOP = ___ , DUST = ___
###
