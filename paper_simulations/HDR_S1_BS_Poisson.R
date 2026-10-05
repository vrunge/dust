
library(changepoint)
library(callr)
library(dust)

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  rates <- rep(c(4, 20), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = rates, type = "poisson")
}

# changepoint::cpt.meanvar(method="BinSeg") can overflow its C stack on
# large/awkward n -- this crashes the R process outright, it is not a
# catchable error (see HDR_chap1_DUST_Poisson.R and
# Section1_introduction_poisson.R, same issue found there for Poisson
# BinSeg/SegNeigh). Unlike the Gaussian file (fpop::multiBinSeg, which
# never crashes), BS is run here in a disposable callr subprocess: if
# that subprocess's C stack overflows, only it dies and callr::r()
# returns NULL to this (healthy) parent, instead of taking down the
# whole 10-rep loop / the rest of the script.
run_BS <- function(data, nb_seg) {
  # A C stack overflow sometimes returns NULL from callr::r() silently, but
  # can also surface as a thrown "subprocess crashed" error in the parent
  # (confirmed empirically -- both failure modes are real, not just one).
  # tryCatch here, not just an is.null() check at the call site, so either
  # one is treated the same way: this n is not usable.
  tryCatch(
    callr::r(function(data, nb_seg) {
      library(changepoint)
      Q <- min(2 * nb_seg + 10, max(nb_seg, floor(length(data) / 3)))
      timing <- system.time(
        res <- changepoint::cpt.meanvar(data, test.stat = "Poisson", method = "BinSeg", penalty = "BIC", Q = Q)
      )
      list(user = timing[["user.self"]], system = timing[["sys.self"]],
           elapsed = timing[["elapsed"]], n_detected = ncpts(res))
    }, args = list(data = data, nb_seg = nb_seg)),
    error = function(e) NULL
  )
}

################################################################################
################################################################################
##### nb_seg = 10
################################################################################
################################################################################
### Re-confirmed this session (crash-boundary bisection, rates 4/8): crash
### boundary is between n=985407 (safe) and n=997500 (crash) -- matches
### the inherited anchor closely. CRASH ceiling, not time ceiling: expect
### an R error well under 10s if the range is pushed too high, not a slow
### run.
################################################################################

nb_seg <- 10

# Store individual runs
results10 <- data.frame()

for (i in seq(from = 9.6, to = 10.0, by = 0.1)) {

  n <- i * 10^5
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    out <- run_BS(data, nb_seg)

    if (is.null(out)) {
      # subprocess crashed (C stack overflow): record as unusable, keep going
      user <- NA; system <- NA; elapsed <- NA; n_detected <- NA
    } else {
      user <- out$user; system <- out$system; elapsed <- out$elapsed; n_detected <- out$n_detected
    }

    results10 <- rbind(
      results10,
      data.frame(
        i = i,
        n = n,
        run = j,
        user = user,
        system = system,
        elapsed = elapsed,
        n_detected = n_detected
      )
    )

    print(j)
  }
}


# Summary statistics
library(dplyr)

summary_results10 <- results10 %>%
  group_by(i, n) %>%
  summarise(
    N_runs = n(),
    N_ok = sum(!is.na(elapsed)),

    mean_elapsed = mean(elapsed, na.rm = TRUE),
    sd_elapsed = sd(elapsed, na.rm = TRUE),
    se_elapsed = sd_elapsed / sqrt(N_ok),

    CI_lower = mean_elapsed -
      qt(0.975, df = N_ok - 1) * se_elapsed,

    CI_upper = mean_elapsed +
      qt(0.975, df = N_ok - 1) * se_elapsed,

    mean_n_detected = mean(n_detected, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(nb_seg = 10)

print(summary_results10)



################################################################################
################################################################################
##### nb_seg = 100
################################################################################
################################################################################
### Narrowed via crash-boundary bisection (rates 4/8): boundary is between
### n=966314 (safe) and n=991115 (crash) -- essentially the SAME boundary
### as nb_seg=10, confirming the crash depends on raw data length n, not
### on nb_seg/Q (Q was 210 here vs 30 at nb_seg=10, 7x different, with no
### effect on the threshold).
################################################################################

nb_seg <- 100

# Store individual runs
results100 <- data.frame()

for (i in seq(from = 9.4, to = 10.0, by = 0.1)) {

  n <- i * 10^5
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    out <- run_BS(data, nb_seg)

    if (is.null(out)) {
      user <- NA; system <- NA; elapsed <- NA; n_detected <- NA
    } else {
      user <- out$user; system <- out$system; elapsed <- out$elapsed; n_detected <- out$n_detected
    }

    results100 <- rbind(
      results100,
      data.frame(
        i = i,
        n = n,
        run = j,
        user = user,
        system = system,
        elapsed = elapsed,
        n_detected = n_detected
      )
    )

    print(j)
  }
}


# Summary statistics
library(dplyr)

summary_results100 <- results100 %>%
  group_by(i, n) %>%
  summarise(
    N_runs = n(),
    N_ok = sum(!is.na(elapsed)),

    mean_elapsed = mean(elapsed, na.rm = TRUE),
    sd_elapsed = sd(elapsed, na.rm = TRUE),
    se_elapsed = sd_elapsed / sqrt(N_ok),

    CI_lower = mean_elapsed -
      qt(0.975, df = N_ok - 1) * se_elapsed,

    CI_upper = mean_elapsed +
      qt(0.975, df = N_ok - 1) * se_elapsed,

    mean_n_detected = mean(n_detected, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(nb_seg = 100)

print(summary_results100)



################################################################################
################################################################################
##### nb_seg = 1000
################################################################################
################################################################################
### Narrowed via crash-boundary bisection (rates 4/8): boundary is between
### n=987164 (safe) and n=1012500 (crash) -- again essentially the same
### n-boundary as nb_seg=10/100 (confirmed across Q=30/210/2010: the crash
### depends on raw data length, not nb_seg/Q). Note elapsed times here are
### notably higher than nb_seg=10/100 at a similar n (up to ~20s) -- BinSeg
### at nb_seg=1000 is close to its own crash ceiling and its time ceiling
### simultaneously, unlike the other two blocks.
################################################################################

nb_seg <- 1000

# Store individual runs
results1000 <- data.frame()

for (i in seq(from = 9.6, to = 10.2, by = 0.1)) {

  n <- i * 10^5
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    out <- run_BS(data, nb_seg)

    if (is.null(out)) {
      user <- NA; system <- NA; elapsed <- NA; n_detected <- NA
    } else {
      user <- out$user; system <- out$system; elapsed <- out$elapsed; n_detected <- out$n_detected
    }

    results1000 <- rbind(
      results1000,
      data.frame(
        i = i,
        n = n,
        run = j,
        user = user,
        system = system,
        elapsed = elapsed,
        n_detected = n_detected
      )
    )

    print(j)
  }
}


# Summary statistics
library(dplyr)

summary_results1000 <- results1000 %>%
  group_by(i, n) %>%
  summarise(
    N_runs = n(),
    N_ok = sum(!is.na(elapsed)),

    mean_elapsed = mean(elapsed, na.rm = TRUE),
    sd_elapsed = sd(elapsed, na.rm = TRUE),
    se_elapsed = sd_elapsed / sqrt(N_ok),

    CI_lower = mean_elapsed -
      qt(0.975, df = N_ok - 1) * se_elapsed,

    CI_upper = mean_elapsed +
      qt(0.975, df = N_ok - 1) * se_elapsed,

    mean_n_detected = mean(n_detected, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(nb_seg = 1000)

print(summary_results1000)



################################################################################
################################################################################
##### Merged summary across all nb_seg
################################################################################
################################################################################

summary_results_BS_Poisson <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_BS_Poisson)

write.csv(summary_results_BS_Poisson,
          file = "paper_simulations/HDR_S1_BS_Poisson_summary.csv",
          row.names = FALSE)
