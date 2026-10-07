# Calibration of the penalty scale factor c0 (beta = 2 c0 log(n)) for the 8
# models of the DUST paper, at n = 1e4:
#   - with no change, DUST must detect no change point;
#   - with k = 1, ..., 9 changes of small amplitude (the two parameter values
#     of Table 2), DUST must detect exactly k change points.
# For each model, every c0 of a grid is tried on the same series.
#
# Self-contained: needs only the package dust.
# Run it from the folder that contains it, in R with source("CALIBRATION.R")
# or in a terminal with: Rscript CALIBRATION.R
# Outputs (same folder): calibration_results.rds (all runs) and
# calibration_summary.csv (one row per model and c0).

library(dust)

## Settings ---------------------------------------------------------------
n           <- 1e4
changes     <- 0:9                                   # true numbers of changes
repetitions <- 100                                   # series per number of changes
c0_grid     <- exp(seq(log(0.01), log(3), length.out = 41))
seed        <- 20261008
size        <- 10                                    # binomial / negative binomial size
# 50 cores, or every core when fewer exist (Windows cannot fork: 1 core).
cores <- if (.Platform$OS.type == "windows") 1 else min(50, parallel::detectCores())

# Table 2 of the paper: the two alternating parameter values and c0.
paper <- list(
  gauss    = list(values = c(0, 1),     c0 = 1),
  poisson  = list(values = c(3, 4),     c0 = 2/3),
  exp      = list(values = c(1, 0.5),   c0 = 3/4),
  geom     = list(values = c(0.5, 0.7), c0 = 2/3),
  bern     = list(values = c(0.5, 0.7), c0 = 2/3),
  binom    = list(values = c(0.5, 0.7), c0 = 1/6),
  negbin   = list(values = c(0.5, 0.7), c0 = 1/10),
  variance = list(values = c(1, 2),     c0 = 1)
)

# k evenly spaced changes, parameters alternating between the two values,
# then data_normalization_1D() as in the paper.
simulate_data <- function(model, k) {
  ends <- round(seq_len(k + 1) * n / (k + 1))
  y <- dataGenerator_1D(chpts = ends, parameters = rep(paper[[model]]$values, length.out = k + 1),
                        nbTrials = size, nbSuccess = size, type = model)
  if (model %in% c("binom", "negbin")) data_normalization_1D(y, type = model, size = size)
  else data_normalization_1D(y, type = model)
}

# Runs task(1), ..., task(n_tasks) on `cores` workers. Each worker gets a new
# task as soon as it finishes the previous one (no idle cores), and the
# progress is printed after every `cores` completed tasks.
run_with_progress <- function(n_tasks, task) {
  start <- Sys.time()
  out <- vector("list", n_tasks)
  report <- function(done) {
    seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    duration <- function(s) if (s < 120) sprintf("%.0f s", s) else sprintf("%.1f min", s / 60)
    cat(sprintf("%s  %d / %d done (%3.0f%%)  elapsed %s, about %s left\n",
                format(Sys.time(), "%H:%M:%S"), done, n_tasks, 100 * done / n_tasks,
                duration(seconds), duration(seconds / done * (n_tasks - done))))
    flush.console()
  }
  if (cores <= 1) {
    for (i in seq_len(n_tasks)) {
      out[[i]] <- try(task(i), silent = TRUE)
      if (i %% 10 == 0 || i == n_tasks) report(i)
    }
    return(out)
  }
  jobs <- list()      # running jobs
  index <- integer()  # task number of each running job
  next_task <- 1
  done <- 0
  while (done < n_tasks) {
    while (length(jobs) < cores && next_task <= n_tasks) {
      i <- next_task
      jobs[[length(jobs) + 1]] <- parallel::mcparallel(try(task(i), silent = TRUE))
      index <- c(index, i)
      next_task <- next_task + 1
    }
    finished <- parallel::mccollect(jobs, wait = FALSE, timeout = 1)
    if (is.null(finished)) next
    pids <- vapply(jobs, function(j) j$pid, numeric(1))
    for (pid in names(finished)) {
      k <- which(pids == as.numeric(pid))
      out[[index[k]]] <- finished[[pid]]
      done <- done + 1
      if (done %% cores == 0 || done == n_tasks) report(done)
    }
    keep <- !(pids %in% as.numeric(names(finished)))
    jobs <- jobs[keep]
    index <- index[keep]
  }
  out
}

## Simulation ---------------------------------------------------------------
# One task = one model, one number of changes, one repetition: one series,
# segmented with every c0 of the grid.
grid <- expand.grid(model = names(paper), k = changes, replicate = seq_len(repetitions),
                    stringsAsFactors = FALSE)
one_task <- function(i) {
  set.seed(seed + i)
  model <- grid$model[i]
  k <- grid$k[i]
  y <- simulate_data(model, k)
  detected <- vapply(c0_grid, function(c0)
    length(dust.1D(y, penalty = 2 * c0 * log(n), model = model)$changepoints) - 1, numeric(1))
  data.frame(model = model, k = k, replicate = grid$replicate[i], c0 = c0_grid, detected = detected)
}
cat(nrow(grid), "series x", length(c0_grid), "penalties on", cores, "cores\n")
runs <- run_with_progress(nrow(grid), one_task)
failed <- vapply(runs, function(r) !is.data.frame(r), logical(1))
if (any(failed)) stop(sum(failed), " tasks failed, for example: ", runs[[which(failed)[1]]])
results <- do.call(rbind, runs)
saveRDS(results, "calibration_results.rds")

## Summary ------------------------------------------------------------------
# For each model and c0:
#   no_change_ok  share of no-change series with no change detected,
#   exact_k       share of series with k = 1..9 changes where exactly k are detected,
#   too_few / too_many  shares of those series with fewer / more than k detected.
summary <- do.call(rbind, lapply(split(results, list(results$model, results$c0), drop = TRUE), function(d) {
  none <- d[d$k == 0, ]
  some <- d[d$k > 0, ]
  data.frame(model = d$model[1], c0 = d$c0[1],
             no_change_ok = mean(none$detected == 0),
             exact_k = mean(some$detected == some$k),
             too_few = mean(some$detected < some$k),
             too_many = mean(some$detected > some$k))
}))
summary <- summary[order(match(summary$model, names(paper)), summary$c0), ]
write.csv(summary, "calibration_summary.csv", row.names = FALSE)

# Recommended c0: among the values with no false detection at all on the
# no-change series, the one detecting exactly k changes most often (the
# smallest such c0 in case of ties).
cat("\nmodel     paper c0  (no change ok, exact k)   calibrated c0  (no change ok, exact k)\n")
for (m in names(paper)) {
  s <- summary[summary$model == m, ]
  at_paper <- s[which.min(abs(log(s$c0) - log(paper[[m]]$c0))), ]
  ok <- s[s$no_change_ok == 1, ]
  best <- if (nrow(ok)) ok[which.max(ok$exact_k), ] else s[which.max(s$no_change_ok), ]
  cat(sprintf("%-9s %6.3f    (%5.1f%%, %5.1f%%)          %6.3f     (%5.1f%%, %5.1f%%)\n",
              m, paper[[m]]$c0, 100 * at_paper$no_change_ok, 100 * at_paper$exact_k,
              best$c0, 100 * best$no_change_ok, 100 * best$exact_k))
}
