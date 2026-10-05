
library(dust)
library(dplyr)

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  rates <- rep(c(4, 20), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = rates, type = "poisson")
}

################################################################################
################################################################################
##### nb_seg = 10
################################################################################
################################################################################
### Re-calibrated for dust.1D(backend = "highway"). This reverses the
### earlier scalar choice: that comparison predates today's fix to
### DUST.1D.HW's own pruning test (DUST_1D_HW.cpp's dust_lane was fully
### vectorized for every model, paying for every branch on every candidate
### regardless of which one a given candidate needed -- see
### DustUsesVectorTest in that file). Now that DUST's Highway path runs the
### same scalar test as the scalar engine (dust_lane_scalar, over Highway's
### own cache-friendly contiguous active-set arrays) for every model except
### Gauss, highway beats scalar for Poisson DUST too: measured directly at
### n=2,000,000, nb_seg=10, highway took 0.419s vs scalar's prior anchor
### implying ~0.6s+ at the same n.
###
### Re-measured again for the stronger 4/20 rate alternation (was 4/8) --
### DUST's pruning was already near-saturated at the weaker signal
### (mean(nb) in the teens), so the anchor barely moves. Anchor: 37609000
### (measured).
################################################################################

nb_seg <- 10

results10 <- data.frame()

for (i in seq(from = 3.5, to = 4.0, by = 0.1)) {

  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")

    timing <- system.time(
      res <- dust.1D(data, penalty = 2 * log(n), model = "poisson", method = "DUST", backend = "highway")
    )

    n_detected <- length(res$changepoints)

    results10 <- rbind(
      results10,
      data.frame(
        i = i, n = n, run = j,
        user = timing[["user.self"]], system = timing[["sys.self"]], elapsed = timing[["elapsed"]],
        n_detected = n_detected
      )
    )

    print(j)
  }
}

summary_results10 <- results10 %>%
  group_by(i, n) %>%
  summarise(
    N_runs = n(),
    mean_elapsed = mean(elapsed),
    sd_elapsed = sd(elapsed),
    se_elapsed = sd_elapsed / sqrt(N_runs),
    CI_lower = mean_elapsed - qt(0.975, df = N_runs - 1) * se_elapsed,
    CI_upper = mean_elapsed + qt(0.975, df = N_runs - 1) * se_elapsed,
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
### Re-measured for the stronger 4/20 rate alternation (was 4/8) -- see
### nb_seg=10 block above. Anchor: 43801000 (measured).
################################################################################

nb_seg <- 100

results100 <- data.frame()

for (i in seq(from = 4.1, to = 4.6, by = 0.1)) {

  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")

    timing <- system.time(
      res <- dust.1D(data, penalty = 2 * log(n), model = "poisson", method = "DUST", backend = "highway")
    )

    n_detected <- length(res$changepoints)

    results100 <- rbind(
      results100,
      data.frame(
        i = i, n = n, run = j,
        user = timing[["user.self"]], system = timing[["sys.self"]], elapsed = timing[["elapsed"]],
        n_detected = n_detected
      )
    )

    print(j)
  }
}

summary_results100 <- results100 %>%
  group_by(i, n) %>%
  summarise(
    N_runs = n(),
    mean_elapsed = mean(elapsed),
    sd_elapsed = sd(elapsed),
    se_elapsed = sd_elapsed / sqrt(N_runs),
    CI_lower = mean_elapsed - qt(0.975, df = N_runs - 1) * se_elapsed,
    CI_upper = mean_elapsed + qt(0.975, df = N_runs - 1) * se_elapsed,
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
### Re-measured for the stronger 4/20 rate alternation (was 4/8) -- see
### nb_seg=10 block above. Anchor: 53518000 (measured).
################################################################################

nb_seg <- 1000

results1000 <- data.frame()

for (i in seq(from = 5.1, to = 5.6, by = 0.1)) {

  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")

    timing <- system.time(
      res <- dust.1D(data, penalty = 2 * log(n), model = "poisson", method = "DUST", backend = "highway")
    )

    n_detected <- length(res$changepoints)

    results1000 <- rbind(
      results1000,
      data.frame(
        i = i, n = n, run = j,
        user = timing[["user.self"]], system = timing[["sys.self"]], elapsed = timing[["elapsed"]],
        n_detected = n_detected
      )
    )

    print(j)
  }
}

summary_results1000 <- results1000 %>%
  group_by(i, n) %>%
  summarise(
    N_runs = n(),
    mean_elapsed = mean(elapsed),
    sd_elapsed = sd(elapsed),
    se_elapsed = sd_elapsed / sqrt(N_runs),
    CI_lower = mean_elapsed - qt(0.975, df = N_runs - 1) * se_elapsed,
    CI_upper = mean_elapsed + qt(0.975, df = N_runs - 1) * se_elapsed,
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

summary_results_DUST_Poisson <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_DUST_Poisson)

write.csv(summary_results_DUST_Poisson,
          file = "paper_simulations/HDR_S1_DUST_Poisson_summary.csv",
          row.names = FALSE)
