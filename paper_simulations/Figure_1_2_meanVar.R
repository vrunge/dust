# Reproduce HDR Figure 1.2 with the original dust meanVar implementation.
# Usage: Rscript Figure_1_2_meanVar.R [replications] [output_directory]
# The original script did not set a seed; this rerun fixes one for reproducibility.

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1L) as.integer(args[1L]) else 10000L
out_dir <- if (length(args) >= 2L) args[2L] else "output/figure_1_2"
library(ggplot2)

# The installed dust package may predate dust.meanVar. Compile exactly the two
# original C++ source files used by dust::dust.meanVar, without editing them.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
dust_src <- normalizePath(file.path(script_dir, "../../dust/src"))
wrapper <- tempfile(fileext = ".cpp")
writeLines(c(
  "#include <Rcpp.h>",
  sprintf('#include "%s"', file.path(dust_src, "2D_Indices.cpp")),
  sprintf('#include "%s"', file.path(dust_src, "2DmeanVar_DUST.cpp")),
  "// [[Rcpp::export]]",
  "Rcpp::List dust_reference(Rcpp::NumericVector data, double penalty, std::string method) {",
  "  std::vector<double> values(data.begin(), data.end());",
  '  DUST_meanVar object(method == "1D" ? "DUST1" : "DUST2",',
  '                      method == "1D" ? "det1" : "det2");',
  "  object.dust(values, Rcpp::Nullable<double>(Rcpp::wrap(penalty)));",
  "  return object.get_partition();",
  "}"
), wrapper)
Rcpp::sourceCpp(wrapper, verbose = FALSE)

stopifnot(length(n_rep) == 1L, !is.na(n_rep), n_rep >= 1L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n <- 10000L
penalty <- 4 * log(n)
seed <- 20261004L
set.seed(seed)

nb1 <- matrix(NA_integer_, nrow = n, ncol = n_rep)
nb2 <- matrix(NA_integer_, nrow = n, ncol = n_rep)
start_time <- proc.time()[[3L]]
for (rep_id in seq_len(n_rep)) {
  y <- rnorm(n)
  nb1[, rep_id] <- dust_reference(y, penalty, "1D")$nb
  nb2[, rep_id] <- dust_reference(y, penalty, "2D")$nb
  if (rep_id %% 100L == 0L || rep_id == n_rep) {
    cat(sprintf("replication %d/%d, elapsed %.1f s\n",
                rep_id, n_rep, proc.time()[[3L]] - start_time))
    flush.console()
  }
}

summarize <- function(x) {
  list(
    mean = rowMeans(x),
    lower = apply(x, 1L, quantile, probs = 0.025, names = FALSE),
    upper = apply(x, 1L, quantile, probs = 0.975, names = FALSE),
    exemplar = x[, 1L]
  )
}
one <- summarize(nb1)
two <- summarize(nb2)
rm(nb1, nb2)
gc(FALSE)

summary <- list(n = n, n_rep = n_rep, penalty = penalty, seed = seed,
                implementation = "original dust source", one = one, two = two)
saveRDS(summary, file.path(out_dir, "figure_1_2_summary.rds"))

pct1 <- 100 * one$mean[n] / n
pct2 <- 100 * two$mean[n] / n
label1 <- sprintf("DUST (1 constr) mean (%.2f%% left at n = %d)", pct1, n)
label2 <- sprintf("DUST (2 constr) mean (%.2f%% left at n = %d)", pct2, n)
time <- seq_len(n)
mean_df <- rbind(
  data.frame(time = time, method = "1D", mean = one$mean,
             lower = one$lower, upper = one$upper),
  data.frame(time = time, method = "2D", mean = two$mean,
             lower = two$lower, upper = two$upper)
)
ex_df <- rbind(
  data.frame(time = time, method = "1D", exemplar = one$exemplar),
  data.frame(time = time, method = "2D", exemplar = two$exemplar)
)
y_min <- min(mean_df$lower, ex_df$exemplar)
y_max <- max(mean_df$upper, ex_df$exemplar)

figure <- ggplot() +
  geom_ribbon(data = mean_df, aes(time, ymin = lower, ymax = upper),
              fill = "grey60", alpha = 0.25) +
  geom_line(data = mean_df, aes(time, mean, colour = method), linewidth = 1) +
  geom_line(data = subset(ex_df, method == "1D"), aes(time, exemplar),
            colour = "#8119FF", linewidth = 0.3, alpha = 0.7) +
  geom_line(data = subset(ex_df, method == "2D"), aes(time, exemplar),
            colour = "#FF19B2", linewidth = 0.3, alpha = 0.7) +
  facet_grid(. ~ method, labeller = as_labeller(c("1D" = "DUST (1 constr)",
                                                "2D" = "DUST (2 constr)"))) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(x = "", y = "number of non-pruned indices", colour = "") +
  scale_colour_manual(values = c("1D" = "#3236FC", "2D" = "#FF334B"),
                      breaks = c("1D", "2D"), labels = c(label1, label2)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "figure_1_2_reproduced.png"), figure,
       width = 12, height = 5.79, dpi = 180, bg = "white")
ggsave(file.path(out_dir, "figure_1_2_reproduced.pdf"), figure,
       width = 12, height = 5.79, device = grDevices::pdf, bg = "white")
cat(sprintf("At n: one constraint %.3f%%, two constraints %.3f%%\n", pct1, pct2))
