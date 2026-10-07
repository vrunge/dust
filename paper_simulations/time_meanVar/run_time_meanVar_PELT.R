# PELT subset of the mean/variance timing experiment from the DUST paper.
# Usage: Rscript run_time_meanVar_PELT.R [repetitions] [output_directory] [workers] [segment_counts] [first_rep]
# segment_counts is a comma-separated subset of 10000,1000,100,10,1.
# The default is 100 repetitions; use a smaller number for a pilot run.

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1L) as.integer(args[1L]) else 100L
if (length(n_rep) != 1L || is.na(n_rep) || n_rep < 1L)
  stop("repetitions must be a positive integer")
first_rep <- if (length(args) >= 5L) as.integer(args[5L]) else 1L
if (length(first_rep) != 1L || is.na(first_rep) || first_rep < 1L)
  stop("first_rep must be a positive integer")

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else ""
script_dir <- if (nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
out_dir <- if (length(args) >= 2L) args[2L] else file.path(script_dir, "pelt")
workers <- if (length(args) >= 3L) as.integer(args[3L]) else
  min(8L, parallel::detectCores(logical = FALSE))
if (length(workers) != 1L || is.na(workers) || workers < 1L)
  stop("workers must be a positive integer")

# Select the PELT build explicitly when it is in a separate library.
benchmark_lib <- Sys.getenv("DUST_BENCHMARK_LIB")
if (nzchar(benchmark_lib)) .libPaths(c(benchmark_lib, .libPaths()))
if (!requireNamespace("dust", quietly = TRUE))
  stop("Install the dust package before running this script")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n <- 1e6L
penalty <- 2 * log(n)
segment_counts <- c(1e4L, 1e3L, 1e2L, 10L, 1L)
selected_segments <- if (length(args) >= 4L)
  as.integer(strsplit(args[4L], ",", fixed = TRUE)[[1L]]) else c(1e4L, 1e3L)
if (!length(selected_segments) || anyNA(selected_segments) ||
    any(!selected_segments %in% segment_counts))
  stop("segment_counts must be a comma-separated subset of 10000,1000,100,10,1")
selected_segments <- unique(selected_segments)
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
  methods <- "PELT"
  if (job$rep %% 2L == 0L) methods <- rev(methods)
  rows <- lapply(methods, function(method) {
    elapsed <- unname(system.time(
      fit <- dust::dust.meanVar(y, penalty = penalty, method = method, backend = "highway")
    )[["elapsed"]])
    data.frame(nb_seg = job$nb_seg, rep = job$rep, seed = job$seed, method = method,
               time_s = elapsed,
               nb_changepoints = max(0L, length(fit$changepoints) - 1L),
               segments_detected = length(fit$changepoints))
  })
  do.call(rbind, rows)
}

jobs <- expand.grid(nb_seg = segment_counts,
                    rep = seq.int(from = first_rep, length.out = n_rep))
jobs$seed <- base_seed + length(segment_counts) * (jobs$rep - 1L) +
  match(jobs$nb_seg, segment_counts)
# Assign seeds before filtering, so each repetition uses the original input.
jobs <- jobs[jobs$nb_seg %in% selected_segments, ]
jobs <- split(jobs, seq_len(nrow(jobs)))
cat(sprintf("Running %d repetitions for each of %d segment counts using %d workers.\n",
            n_rep, length(selected_segments), workers))
start <- proc.time()[["elapsed"]]
result <- vector("list", length(jobs))
for (first in seq.int(1L, length(jobs), by = workers)) {
  batch <- first:min(first + workers - 1L, length(jobs))
  result[batch] <- parallel::mclapply(jobs[batch], one_run, mc.cores = workers,
                                    mc.preschedule = FALSE)
  if (any(vapply(result[batch], inherits, logical(1L), "try-error")))
    stop("A PELT simulation failed; inspect the worker error above")
  write.csv(do.call(rbind, result), file.path(out_dir, "timings.csv"), row.names = FALSE)
  cat(sprintf("Completed %d/%d fits in %.1f seconds.\n", max(batch), length(jobs),
              proc.time()[["elapsed"]] - start))
}

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
saveRDS(list(n = n, penalty = penalty, repetitions = n_rep, first_rep = first_rep,
             base_seed = base_seed, workers = workers, backend = "highway",
             library = find.package("dust"), timings = timings, summary = summary,
             elapsed_s = proc.time()[["elapsed"]] - start,
             dust_version = as.character(utils::packageVersion("dust"))),
        file.path(out_dir, "timings.rds"))

cat(sprintf("Completed in %.1f seconds. PELT results saved in %s\n",
            proc.time()[["elapsed"]] - start, normalizePath(out_dir)))
print(summary, row.names = FALSE)
