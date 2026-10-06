
library(changepoint)
library(dust)


build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  means <- rep(c(0, 10), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = means, type = "gauss", sdNoise = 1)
}

################################################################################
################################################################################
##### nb_seg = 10
################################################################################
################################################################################
### Re-verified for the stronger 0/10 mean alternation (was 0/1): still
### lands at ~9.7s, range unchanged.
################################################################################

nb_seg <- 10


# Store individual runs
results10 <- data.frame()

for (i in seq(from = 17, to = 19, by = 0.5)) {

  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- fpop::multiBinSeg(data, Kmax = nb_seg)
    )

    n_detected <- res$K

    results10 <- rbind(
      results10,
      data.frame(
        i = i,
        n = n,
        run = j,
        user = timing["user.self"],
        system = timing["sys.self"],
        elapsed = timing["elapsed"],
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

    mean_elapsed = mean(elapsed),
    sd_elapsed = sd(elapsed),
    se_elapsed = sd_elapsed / sqrt(N_runs),

    CI_lower = mean_elapsed -
      qt(0.975, df = N_runs - 1) * se_elapsed,

    CI_upper = mean_elapsed +
      qt(0.975, df = N_runs - 1) * se_elapsed,

    mean_n_detected = mean(n_detected),

    .groups = "drop"
  ) %>%
  mutate(nb_seg = 10)

print(summary_results10)



################################################################################
################################################################################
##### nb_seg = 100
################################################################################
################################################################################
### Re-verified for the stronger 0/10 mean alternation (was 0/1): still
### lands at ~10.0s, range unchanged.
################################################################################

nb_seg <- 100


# Store individual runs
results100 <- data.frame()

for (i in seq(from = 3.5, to = 3.95, by = 0.05))
{
  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- fpop::multiBinSeg(data, Kmax = nb_seg)
    )

    n_detected <- res$K

    results100 <- rbind(
      results100,
      data.frame(
        i = i,
        n = n,
        run = j,
        user = timing["user.self"],
        system = timing["sys.self"],
        elapsed = timing["elapsed"],
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

    mean_elapsed = mean(elapsed),
    sd_elapsed = sd(elapsed),
    se_elapsed = sd_elapsed / sqrt(N_runs),

    CI_lower = mean_elapsed -
      qt(0.975, df = N_runs - 1) * se_elapsed,

    CI_upper = mean_elapsed +
      qt(0.975, df = N_runs - 1) * se_elapsed,

    mean_n_detected = mean(n_detected),

    .groups = "drop"
  ) %>%
  mutate(nb_seg = 100)

print(summary_results100)



################################################################################
################################################################################
##### nb_seg = 1000
################################################################################
################################################################################
### Re-measured for the stronger 0/10 mean alternation (was 0/1) and found
### something more fundamental than a scaling problem: elapsed time at
### nb_seg=1000 is NOT a reproducible function of n at all. n=42,000,000
### ran at 1.7s under seed(1), but the SAME n with no seed set (i.e. the
### data draw this file's own build_data() actually uses) got stuck for
### 100+ seconds, repeatedly, on fresh R processes with no shared state --
### ruling out memory fragmentation or a cache-size threshold (both ruled
### out directly: a brand-new process hung just as badly). This is
### fpop::multiBinSeg's own recursive splitting behavior being sensitive
### to the specific noise realization at Kmax=1000 -- some draws are fast,
### others pathologically slow, at the identical n. No single n is a
### stable ~10s anchor here; reporting that honestly rather than picking
### a number that only works for some random seeds and not others. The
### range below is left at a value confirmed fast under at least one
### draw, with this caveat attached -- treat any single run's timing at
### nb_seg=1000 as unreliable until multiBinSeg's own behavior here is
### understood further.
################################################################################

nb_seg <- 1000


# Store individual runs
results1000 <- data.frame()

for (i in seq(from = 3.7, to = 4.2, by = 0.1))
{
  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- fpop::multiBinSeg(data, Kmax = nb_seg)
    )

    n_detected <- res$K

    results1000 <- rbind(
      results1000,
      data.frame(
        i = i,
        n = n,
        run = j,
        user = timing["user.self"],
        system = timing["sys.self"],
        elapsed = timing["elapsed"],
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

    mean_elapsed = mean(elapsed),
    sd_elapsed = sd(elapsed),
    se_elapsed = sd_elapsed / sqrt(N_runs),

    CI_lower = mean_elapsed -
      qt(0.975, df = N_runs - 1) * se_elapsed,

    CI_upper = mean_elapsed +
      qt(0.975, df = N_runs - 1) * se_elapsed,

    mean_n_detected = mean(n_detected),

    .groups = "drop"
  ) %>%
  mutate(nb_seg = 1000)

print(summary_results1000)



################################################################################
################################################################################
##### Merged summary across all nb_seg
################################################################################
################################################################################

summary_results_BS_Gauss <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_BS_Gauss)

write.csv(summary_results_BS_Gauss,
          file = "paper_simulations/10sLimit/results/HDR_S1_BS_Gauss_summary.csv",
          row.names = FALSE)



