
###
### Section 1. INTRODUCTION (Gaussian model)
###

#########
# system time, limit = 10s
# Three experiments: 10, 100 and 1000 equal-length segments (9, 99, 999
# true changepoints), means alternating 0/1. For each, find the largest
# n every algorithm can segment within a 10s budget.
#########
# BS
# OP
# PELT
# FPOP
# DUST (dust.1D(..., backend = "highway"))

# install.packages("fpop", repos="http://R-Forge.R-project.org")
library(fpop)
library(dust)

# Re-run on Apple M4 Pro (14-core: 10 performance + 4 efficiency), 24GB,
# macOS 26.6.2 -- NOT the MacBook Pro M1 described in the paper's current
# footnote (8-core: 4+4, 16GB, macOS Tahoe 26.1). Update the footnote's
# hardware description before using these numbers in the paper.
#
# OP and PELT go through dust with backend = "scalar". OP is the
# unpruned dynamic-programming reference; PELT applies its pruning rule.
#
# BS and FPOP need fpop::multiBinSeg / fpop::Fpop: dust doesn't
# implement either (not DUST-family pruning algorithms). BS's Kmax is
# set to (nb_seg - 1) in every experiment, matching the true number of
# changepoints, as in the original 10-segment run.
#
# DUST uses dust.1D(..., backend = "highway"), with scalar fallback if
# Highway is unavailable in the installed package. Bisection search replaces a manual
# scan (as the original 10-segment experiment used): automatically
# brackets and refines the max-n threshold, so the same code handles
# all three (very different) regimes without hand-tuned search ranges
# per experiment.

################################################################################
##### Generic bisection search for "max n within budget seconds"
################################################################################

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  means <- rep(c(0, 1), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = means, type = "gauss", sdNoise = 1)
}

time_it <- function(run_fun, n, nb_seg) {
  data <- build_data(n, nb_seg)
  system.time(run_fun(data, n))[["elapsed"]]
}

find_max_n <- function(run_fun, nb_seg, n_lo, n_hi, label, budget = 10, tol = 0.02) {
  t_lo <- time_it(run_fun, n_lo, nb_seg)
  while (t_lo > budget) { n_hi <- n_lo; n_lo <- n_lo / 2; t_lo <- time_it(run_fun, n_lo, nb_seg) }
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

run_BS   <- function(data, n, nb_seg) fpop::multiBinSeg(data, Kmax = nb_seg - 1)
run_OP   <- function(data, n) dust::dust.1D(data, penalty = 2*log(n), model = "gauss", method = "OP", backend = "scalar")
run_PELT <- function(data, n) dust::dust.1D(data, penalty = 2*log(n), model = "gauss", method = "PELT", backend = "scalar")
run_FPOP <- function(data, n) fpop::Fpop(x = data, lambda = 2*log(n))
run_DUST <- function(data, n) dust::dust.1D(data, 2*log(n), model = "gauss", method = "DUST", backend = "highway")

results <- data.frame(nb_seg = integer(0), algo = character(0), max_n = numeric(0))

for (nb_seg in c(10, 100, 1000)) {
  n_BS   <- find_max_n(function(d,n) run_BS(d,n,nb_seg),   nb_seg, n_lo = 1e6, n_hi = 2e8, label = "BS")
  n_OP   <- find_max_n(run_OP,   nb_seg, n_lo = 1e4, n_hi = 2e6,  label = "OP")
  n_PELT <- find_max_n(run_PELT, nb_seg, n_lo = 1e4, n_hi = 2e6,  label = "PELT")
  n_FPOP <- find_max_n(run_FPOP, nb_seg, n_lo = 1e6, n_hi = 1e8,  label = "FPOP")
  n_DUST <- find_max_n(run_DUST, nb_seg, n_lo = 1e6, n_hi = 1.5e8, label = "DUST")

  results <- rbind(results, data.frame(
    nb_seg = nb_seg,
    algo   = c("BS", "OP", "PELT", "FPOP", "DUST"),
    max_n  = c(n_BS, n_OP, n_PELT, n_FPOP, n_DUST)
  ))
}

################################################################################
cat("================ SUMMARY (max n within 10s) ================\n")
print(results, row.names = FALSE)

###
### nb_seg =   10 : BS = 172015625 , OP =    285957 , PELT =    285957 , FPOP = 56687500 , DUST =  83648438
### nb_seg =  100 : BS =  37146484 , OP =    911719 , PELT =    911719 , FPOP = 62101563 , DUST =  92960938
### nb_seg = 1000 : BS =   4012207 , OP =   2343750 , PELT =   2750000 , FPOP = 72156250 , DUST = 106929688
###
### BS collapses as segments grow (O(Kmax*n): more splits to find costs
### more full scans), while OP/PELT/FPOP/DUST all *increase* -- more
### frequent genuine changes prune more aggressively (smaller active
### sets), the same pattern seen throughout this session at smaller n.
### At 1000 segments DUST (107M) overtakes BS (4M) by ~25x, flipping the
### "BS is the fast approximate baseline" framing from the 10-segment
### case entirely.
###
