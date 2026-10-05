
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
### Re-measured for the stronger 4/20 rate alternation (was 4/8). PELT's
### active set grew substantially under the stronger signal (mean(nb) now
### in the thousands), yet elapsed time landed close to the same ~10s
### target anyway -- see HDR_S1_PELT_dust_Gauss.R for the same effect.
### Anchor: 208000 (measured).
################################################################################

nb_seg <- 10

results10 <- data.frame()

for (i in seq(from = 1.95, to = 2.2, by = 0.05)) {

  n <- i * 10^5
  print(n)

  for (j in 1:10) {

    data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")

    timing <- system.time(
      res <- dust.1D(data, penalty = 2 * log(n), model = "poisson", method = "PELT", backend = "scalar")
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
### nb_seg=10 block above. Anchor: 646251 (measured).
################################################################################

nb_seg <- 100

results100 <- data.frame()

for (i in seq(from = 6.05, to = 6.80, by = 0.15)) {

  n <- i * 10^5
  print(n)

  for (j in 1:10) {

    data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")

    timing <- system.time(
      res <- dust.1D(data, penalty = 2 * log(n), model = "poisson", method = "PELT", backend = "scalar")
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
### nb_seg=10 block above. Anchor: 1930800 (measured).
################################################################################

nb_seg <- 1000

results1000 <- data.frame()

for (i in seq(from = 1.8, to = 2.05, by = 0.05)) {

  n <- i * 10^6
  print(n)

  for (j in 1:10) {

    data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")

    timing <- system.time(
      res <- dust.1D(data, penalty = 2 * log(n), model = "poisson", method = "PELT", backend = "scalar")
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

summary_results_PELT_dust_Poisson <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_PELT_dust_Poisson)

write.csv(summary_results_PELT_dust_Poisson,
          file = "paper_simulations/HDR_S1_PELT_dust_Poisson_summary.csv",
          row.names = FALSE)
