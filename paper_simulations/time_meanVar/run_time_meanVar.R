# Reproduce the mean/variance timing experiment from the DUST paper.
# Usage: Rscript run_time_meanVar.R [repetitions] [output_directory] [workers]
# The default is 100 repetitions; use a smaller number for a pilot run.

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1L) as.integer(args[1L]) else 100L
if (length(n_rep) != 1L || is.na(n_rep) || n_rep < 1L)
  stop("repetitions must be a positive integer")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
out_dir <- if (length(args) >= 2L) args[2L] else script_dir
workers <- if (length(args) >= 3L) as.integer(args[3L]) else
  min(8L, parallel::detectCores(logical = FALSE))
if (length(workers) != 1L || is.na(workers) || workers < 1L)
  stop("workers must be a positive integer")

if (!requireNamespace("dust", quietly = TRUE))
  stop("Install the dust package before running this script")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n <- 1e6L
penalty <- 2 * log(n)
segment_counts <- c(1e4L, 1e3L, 1e2L, 10L, 1L)
base_seed <- 123L

# This is the generator used by the original timing script in dustOld. At
# every segment boundary it changes either the mean or the variance.
make_series <- function(nb_seg, seed) {
  set.seed(seed)
  segment_length <- n / nb_seg
  means <- numeric(nb_seg)
  variances <- numeric(nb_seg)
  means[1L] <- 4
  variances[1L] <- 1
  if (nb_seg > 1L) {
    for (k in 2:nb_seg) {
      means[k] <- means[k - 1L]
      variances[k] <- variances[k - 1L]
      if (runif(1L) < 0.5) {
        means[k] <- -means[k]
      } else {
        variances[k] <- if (variances[k] == 1) 16 else 1
      }
    }
  }
  y <- numeric(n)
  for (k in seq_len(nb_seg)) {
    first <- as.integer((k - 1L) * segment_length) + 1L
    last <- as.integer(k * segment_length)
    y[first:last] <- rnorm(last - first + 1L, means[k], sqrt(variances[k]))
  }
  y
}

one_run <- function(job) {
  y <- make_series(job$nb_seg, job$seed)
  methods <- c("1D", "2D")
  if (job$rep %% 2L == 0L) methods <- rev(methods)
  rows <- lapply(methods, function(method) {
    elapsed <- unname(system.time(
      fit <- dust::dust.meanVar(y, penalty = penalty, method = method)
    )[["elapsed"]])
    data.frame(nb_seg = job$nb_seg, rep = job$rep, method = method,
               time_s = elapsed,
               nb_changepoints = max(0L, length(fit$changepoints) - 1L),
               segments_detected = length(fit$changepoints))
  })
  do.call(rbind, rows)
}

jobs <- expand.grid(nb_seg = segment_counts, rep = seq_len(n_rep))
jobs$seed <- base_seed + seq_len(nrow(jobs))
jobs <- split(jobs, seq_len(nrow(jobs)))
cat(sprintf("Running %d repetitions for each of %d segment counts using %d workers.\n",
            n_rep, length(segment_counts), workers))
start <- proc.time()[["elapsed"]]
result <- parallel::mclapply(jobs, one_run, mc.cores = workers,
                             mc.preschedule = FALSE)

timings <- do.call(rbind, result)
summary <- aggregate(time_s ~ nb_seg + method, timings,
                     function(x) c(mean = mean(x), sd = sd(x)))
summary <- data.frame(nb_seg = summary$nb_seg, method = summary$method,
                      mean_s = summary$time_s[, "mean"],
                      sd_s = summary$time_s[, "sd"])
summary$mean_segments <- vapply(seq_len(nrow(summary)), function(i) {
  mean(timings$segments_detected[timings$nb_seg == summary$nb_seg[i] &
                                  timings$method == summary$method[i]])
}, numeric(1L))

write.csv(timings, file.path(out_dir, "timings.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
saveRDS(list(n = n, penalty = penalty, repetitions = n_rep,
             base_seed = base_seed, timings = timings, summary = summary,
             elapsed_s = proc.time()[["elapsed"]] - start,
             dust_version = as.character(utils::packageVersion("dust"))),
        file.path(out_dir, "timings.rds"))

source("write_table.R", local = TRUE)
write_meanvar_table(out_dir)
cat(sprintf("Completed in %.1f seconds. Results and standalone table saved in %s\n",
            proc.time()[["elapsed"]] - start, normalizePath(out_dir)))



