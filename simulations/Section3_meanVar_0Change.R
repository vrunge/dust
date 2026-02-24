
################################################################################
# simple estimation of the number of non-pruned indices for n = 10^6

nbIndices <- NULL
for(i in 1:10)
{
  print(i)
  print(system.time(res <- dust.meanVar(rnorm(10^6), penalty = 4*log(10^6), method = "det_DUST1")))
  nbIndices <- c(nbIndices, res$nb[10^6])
}

nbIndices
mean(nbIndices)/10^6 * 100
sd(nbIndices)/10^6 * 100


################################################################################


#devtools::install_github("vrunge/dust")
library(dust)
library(ggplot2)


# ---- One simulation returning nb for BOTH methods ----
oneSimu_nb <- function(n, pen = 4*log(n))
{
  y <- dataGenerator_meanVar(
    chpts = (1:50) * (n/50),
    means = rep(0, 50),  # your sample(c(0,0),...) is always 0
    sds   = rep(1, 50)   # your sample(c(1,1),...) is always 1
  )

  res1 <- dust.meanVar(y, penalty = pen, method = "det_DUST1")
  res2 <- dust.meanVar(y, penalty = pen, method = "det2_DUST2")

  # Progress marker (like your cat("0"))
  cat("0")

  list(nb1 = res1$nb, nb2 = res2$nb)
}

# ---- Parameters ----
n_steps <- 1e4
n_rep   <- 10000
pen     <- 4 * log(n_steps)

# ---- Run replications ----
sim_list <- replicate(n_rep, oneSimu_nb(n_steps, pen), simplify = FALSE)
cat("\n")  # newline after progress zeros

# Convert list -> matrices (n_steps x n_rep)
nb1_mat <- do.call(cbind, lapply(sim_list, `[[`, "nb1"))
nb2_mat <- do.call(cbind, lapply(sim_list, `[[`, "nb2"))

##############################################################################
##############################################################################

# Summary per time t
mean1  <- rowMeans(nb1_mat)
lower1 <- apply(nb1_mat, 1, quantile, probs = 0.025, names = FALSE)
upper1 <- apply(nb1_mat, 1, quantile, probs = 0.975, names = FALSE)

mean2  <- rowMeans(nb2_mat)
lower2 <- apply(nb2_mat, 1, quantile, probs = 0.025, names = FALSE)
upper2 <- apply(nb2_mat, 1, quantile, probs = 0.975, names = FALSE)

time <- seq_len(n_steps)

df_mean <- data.frame(
  time = time,
  mean1 = mean1, low1 = lower1, up1 = upper1,
  mean2 = mean2, low2 = lower2, up2 = upper2
)

# Optional exemplar curves: take replication 1 for each method
df_ex <- data.frame(
  time = time,
  ex1 = nb1_mat[, 1],
  ex2 = nb2_mat[, 1]
)

y_min <- min(df_mean$low1, df_mean$low2, df_ex$ex1, df_ex$ex2, na.rm = TRUE)
y_max <- max(df_mean$up1,  df_mean$up2,  df_ex$ex1, df_ex$ex2, na.rm = TRUE)


##############################################################################
##############################################################################


library(ggplot2)

# --- % left at time n FROM THE MEAN CURVE ---
pct1 <- 100 * (df_mean$mean1[n_steps] / n_steps)
pct2 <- 100 * (df_mean$mean2[n_steps] / n_steps)

lab1 <- sprintf("DUST (1 constr) mean (%.2f%% left at n = %.0f)", pct1, n_steps)
lab2 <- sprintf("DUST (2 constr) mean (%.2f%% left at n = %.0f)", pct2, n_steps)

# ---- long df for mean + CI (for facets) ----
df_long_mean <- rbind(
  data.frame(time = df_mean$time,
             method = "det_DUST1",
             method_leg = lab1,
             mean = df_mean$mean1, lower = df_mean$low1, upper = df_mean$up1),
  data.frame(time = df_mean$time,
             method = "det2_DUST2",
             method_leg = lab2,
             mean = df_mean$mean2, lower = df_mean$low2, upper = df_mean$up2)
)

# ---- long df for exemplars (no legend, grey) ----
df_long_ex1 <- rbind(
  data.frame(time = df_ex$time, method = "det_DUST1", ex = df_ex$ex1)
)
df_long_ex2 <- rbind(
  data.frame(time = df_ex$time, method = "det2_DUST2", ex = df_ex$ex2)
)



ggplot() +
  geom_ribbon(data = df_long_mean,
              aes(x = time, ymin = lower, ymax = upper),
              fill = "grey60", alpha = 0.25) +
  geom_line(data = df_long_mean,
            aes(x = time, y = mean, color = method),
            linewidth = 1) +
  geom_line(data = df_long_ex1,
            aes(x = time, y = ex),
            color = "#8119FF", linewidth = 0.3, alpha = 0.7,
            show.legend = FALSE) +
  geom_line(data = df_long_ex2,
            aes(x = time, y = ex),
            color = "#FF19B2", linewidth = 0.3, alpha = 0.7,
            show.legend = FALSE) +
  facet_grid(. ~ method,
             labeller = as_labeller(c(det_DUST1="DUST (1 constr)",
                                      det2_DUST2="DUST (2 constr)"))) +
  coord_cartesian(ylim = c(y_min, y_max)) +
  labs(x="", y="number of non-pruned indices", color="") +
  scale_color_manual(values = c("det_DUST1"="#3236FC","det2_DUST2"="#FF334B"),
                     breaks = c("det_DUST1","det2_DUST2"),
                     labels = c("det_DUST1"=lab1,"det2_DUST2"=lab2)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

10

(10^4*10^4*(10^4+1)/2)/sum(nb1_mat)
(10^4*10^4*(10^4+1)/2)/sum(nb2_mat)


