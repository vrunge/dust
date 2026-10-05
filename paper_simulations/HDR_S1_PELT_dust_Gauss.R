
library(dust)
library(dplyr)

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
### Re-calibrated for dust.1D(backend = "highway") (faster than the scalar backend for Gaussian
### PELT, confirmed this session).
###
### Re-measured for the stronger 0/10 mean alternation (was 0/1). Unlike
### DUST, PELT's active set grew substantially under the stronger signal
### (mean(nb) now in the thousands, up from the tens) -- yet elapsed time
### landed close to the same ~10s target anyway, not something chased
### further here (the job is finding n for ~10s, not explaining why).
### Anchor: 377386 (measured).
################################################################################

nb_seg <- 10

results10 <- data.frame()

for (i in seq(from = 3.5, to = 4.0, by = 0.1)) {

  n <- i * 10^5
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- dust.1D(data, 2 * log(n), model = "gauss", method = "PELT", backend = "highway")
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
### Re-measured for the stronger 0/10 mean alternation (was 0/1) -- see
### nb_seg=10 block above. Anchor: 1182829 (measured).
################################################################################

nb_seg <- 100

results100 <- data.frame()

for (i in seq(from = 1.08, to = 1.28, by = 0.04)) {

  n <- i * 10^6
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- dust.1D(data, 2 * log(n), model = "gauss", method = "PELT", backend = "highway")
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
### Re-measured for the stronger 0/10 mean alternation (was 0/1) -- see
### nb_seg=10 block above. Anchor: 3707000 (measured).
################################################################################

nb_seg <- 1000

results1000 <- data.frame()

for (i in seq(from = 3.5, to = 4.0, by = 0.1)) {

  n <- i * 10^6
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- dust.1D(data, 2 * log(n), model = "gauss", method = "PELT", backend = "highway")
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

summary_results_PELT_dust_Gauss <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_PELT_dust_Gauss)

write.csv(summary_results_PELT_dust_Gauss,
          file = "paper_simulations/HDR_S1_PELT_dust_Gauss_summary.csv",
          row.names = FALSE)
