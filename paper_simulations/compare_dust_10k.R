# Reproducible small-data comparison: 8 models x {DUST, DUSTib} x
# {scalar, Highway}. The same seeded data are used by all four runs per model.
# Run with an installed dust package:
#   Rscript paper_simulations/compare_dust_10k.R

library(dust)

n <- 10000L
chpts <- c(2500L, 5000L, 7500L, n)
base_penalty <- 2 * log(n)
known_size <- 10L
models <- c("gauss", "exp", "poisson", "geom", "bern", "binom",
            "negbin", "variance")

make_data <- function(model, index) {
  set.seed(20261004L + index)
  switch(model,
    gauss = dataGenerator_1D(chpts, c(0, 2, -1, 1), sdNoise = 1,
                             type = "gauss"),
    exp = dataGenerator_1D(chpts, c(2, 0.5, 3, 1), type = "exp"),
    poisson = dataGenerator_1D(chpts, c(2, 8, 3, 6), type = "poisson"),
    geom = dataGenerator_1D(chpts, c(0.8, 0.2, 0.7, 0.3), type = "geom"),
    bern = dataGenerator_1D(chpts, c(0.2, 0.8, 0.3, 0.7), type = "bern"),
    binom = dataGenerator_1D(chpts, c(0.2, 0.8, 0.3, 0.7),
                             nbTrials = known_size, type = "binom") / known_size,
    negbin = dataGenerator_1D(chpts, c(0.7, 0.2, 0.6, 0.3),
                              nbSuccess = known_size, type = "negbin") / known_size,
    variance = dataGenerator_1D(chpts, c(1, 3, 0.8, 2), type = "variance")
  )
}

as_index <- function(x) as.integer(unlist(x, use.names = FALSE))
max_abs <- function(x) if (length(x)) max(abs(x)) else 0

# A warm-up call is excluded. Each reported time is a median over five timed
# batches; batch length is chosen from a short pilot to make millisecond-scale
# calls measurable. Timing repeats reuse the same observations.
median_ms <- function(y, penalty, model, method, backend) {
  run <- function() dust.1D(y, penalty, model, method, backend = backend)
  invisible(run())
  pilot <- system.time(invisible(run()))[["elapsed"]]
  repetitions <- min(100L, max(3L, ceiling(0.15 / max(pilot, 0.001))))
  batches <- numeric(5)
  for (j in seq_along(batches)) {
    gc(verbose = FALSE)
    batches[j] <- system.time(
      for (k in seq_len(repetitions)) invisible(run())
    )[["elapsed"]] * 1000 / repetitions
  }
  median(batches)
}

rows <- list()
results <- list()
for (i in seq_along(models)) {
  model <- models[i]
  y <- make_data(model, i)
  stopifnot(length(y) == n, all(is.finite(y)))
  penalty <- base_penalty / if (model %in% c("binom", "negbin")) known_size else 1

  for (method in c("DUST", "DUSTib")) {
    for (backend in c("scalar", "highway")) {
      key <- paste(model, method, backend, sep = "/")
      result <- dust.1D(y, penalty, model, method, backend = backend)
      results[[key]] <- result
      rows[[length(rows) + 1L]] <- data.frame(
        model = model, method = method, backend = backend, n = n,
        penalty = penalty, median_ms = median_ms(y, penalty, model, method, backend),
        changepoints = length(as_index(result$changepoints)),
        active_last = length(as_index(result$lastIndexSet)),
        mean_nb = mean(as_index(result$nb)),
        final_costQ = tail(as.numeric(result$costQ), 1)
      )
      cat(key, "complete\n")
    }
  }
}

summary <- do.call(rbind, rows)
compare <- list()
for (model in models) {
  for (method in c("DUST", "DUSTib")) {
    scalar <- results[[paste(model, method, "scalar", sep = "/")]]
    highway <- results[[paste(model, method, "highway", sep = "/")]]
    cost_diff <- as.numeric(highway$costQ) - as.numeric(scalar$costQ)
    nb_diff <- as_index(highway$nb) - as_index(scalar$nb)
    compare[[length(compare) + 1L]] <- data.frame(
      model = model, method = method,
      changepoints_equal = identical(as_index(highway$changepoints),
                                     as_index(scalar$changepoints)),
      costQ_max_abs_diff = max_abs(cost_diff),
      costQ_max_rel_diff = max_abs(cost_diff / (1 + abs(as.numeric(scalar$costQ)))),
      lastIndexSet_equal = identical(as_index(highway$lastIndexSet),
                                    as_index(scalar$lastIndexSet)),
      nb_equal = identical(as_index(highway$nb), as_index(scalar$nb)),
      nb_max_abs_diff = max_abs(nb_diff)
    )
  }
}
backend_comparison <- do.call(rbind, compare)

method_comparison <- list()
for (model in models) for (backend in c("scalar", "highway")) {
  dust <- results[[paste(model, "DUST", backend, sep = "/")]]
  dustib <- results[[paste(model, "DUSTib", backend, sep = "/")]]
  method_comparison[[length(method_comparison) + 1L]] <- data.frame(
    model = model, backend = backend,
    changepoints_equal = identical(as_index(dust$changepoints), as_index(dustib$changepoints)),
    costQ_max_abs_diff = max_abs(as.numeric(dust$costQ) - as.numeric(dustib$costQ)),
    lastIndexSet_equal = identical(as_index(dust$lastIndexSet), as_index(dustib$lastIndexSet)),
    nb_equal = identical(as_index(dust$nb), as_index(dustib$nb))
  )
}
method_comparison <- do.call(rbind, method_comparison)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
output_dir <- if (length(script_arg))
  dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else "paper_simulations"
if (!dir.exists(output_dir)) output_dir <- "."
write.csv(summary, file.path(output_dir, "compare_dust_10k_results.csv"), row.names = FALSE)
write.csv(backend_comparison,
          file.path(output_dir, "compare_dust_10k_backend_agreement.csv"), row.names = FALSE)
write.csv(method_comparison,
          file.path(output_dir, "compare_dust_10k_method_agreement.csv"), row.names = FALSE)

print(summary, row.names = FALSE)
print(backend_comparison, row.names = FALSE)
print(method_comparison, row.names = FALSE)
