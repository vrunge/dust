
library(fpop)
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
### Range from HDR_chap1_DUST_Gauss.R -- anchor: 56687500 (measured).
### Re-verified for the stronger 0/10 mean alternation (was 0/1): still
### lands at ~10.1s, anchor unchanged.
################################################################################

nb_seg <- 10

results10 <- data.frame()

for (i in seq(from = 5.4, to = 5.95, by = 0.1)) {

  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- fpop::Fpop(x = data, lambda = 2 * log(n))
    )

    n_detected <- length(res$t.est)

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
### Range from HDR_chap1_DUST_Gauss.R -- anchor: 62101563 (measured).
### Re-verified for the stronger 0/10 mean alternation (was 0/1): still
### lands at ~10.0s, anchor unchanged.
################################################################################

nb_seg <- 100

results100 <- data.frame()

for (i in seq(from = 5.9, to = 6.5, by = 0.1)) {

  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- fpop::Fpop(x = data, lambda = 2 * log(n))
    )

    n_detected <- length(res$t.est)

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
### Re-measured for the stronger 0/10 mean alternation (was 0/1) -- old
### anchor (72156250) now ran at 11.47s. Anchor: 67369000 (measured).
################################################################################

nb_seg <- 1000

results1000 <- data.frame()

for (i in seq(from = 6.4, to = 7.0, by = 0.1)) {

  n <- i * 10^7
  print(n)

  for (j in 1:10) {

    data <- build_data(n, nb_seg)

    timing <- system.time(
      res <- fpop::Fpop(x = data, lambda = 2 * log(n))
    )

    n_detected <- length(res$t.est)

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

summary_results_FPOP_Gauss <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_FPOP_Gauss)

write.csv(summary_results_FPOP_Gauss,
          file = "paper_simulations/10sLimit/results/HDR_S1_FPOP_Gauss_summary.csv",
          row.names = FALSE)
