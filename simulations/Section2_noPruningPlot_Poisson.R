library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------
# Parameters + root solve
# -----------------------------
n <- 20
Y <- 30

beta <- Y * log(n/(n-1)) * 9.95/10 * n

poisson_zero_2 <- function(t, n, x, Y, beta) {
  (t/n) * x * log(x) - Y * log(Y) + (Y - (t/n) * x) * log((Y - (t/n) * x) / (1 - (t/n))) - beta/n
}

y <- numeric(n)
for (t in 1:(n - 1)) {
  f_t <- function(x) poisson_zero_2(t, n, x, Y, beta)
  y[t] <- uniroot(f_t, interval = c(0.0001, Y))$root
}
y[n] <- Y

resPoisson <- c(y[1], diff(seq_len(n) * y))
mean(resPoisson)

# -----------------------------
# Plot 1: resPoisson with ggplot2
# -----------------------------
df_res <- tibble(t = seq_len(n), value = resPoisson)

p1 <- ggplot(df_res, aes(x = t, y = value)) +
  geom_line(linewidth = 1.2, color = "red") +
  geom_point(size = 2.6, color = "red") +
  scale_x_continuous(breaks = seq_len(n)) +
  coord_cartesian(ylim = range(resPoisson)) +
  labs(x = "t", y = "resPoisson") +
  theme_minimal(base_size = 14)

print(p1)

# -----------------------------
# Rounding / cumy (as in your code)
# -----------------------------
round(resPoisson)
ceiling(resPoisson)

cumy <- round(cumsum(resPoisson))

# -----------------------------
# Plot 2: Cost curves with ggplot2
# -----------------------------
l <- 0.3
r <- 3
a <- 3.3
b <- 4.5

entropy_term <- function(u) ifelse(u == 0, 0, u - u * log(u))

# mini (for y-limits), same logic as your function
u1 <- cumy[1]
const1 <- entropy_term(u1)
x0 <- log((cumy[n] - cumy[1]) / (n - 1))
mini <- (n - 1) * exp(x0) - (cumy[n] - cumy[1]) * x0 + const1

x_grid <- seq(a, b, length.out = 1000)

# Build a tidy data frame of all curves
cost_df <- bind_rows(
  lapply(1:(n - 1), function(t) {
    u <- cumy[t] / t
    const <- t * entropy_term(u)
    tibble(
      x = x_grid,
      curve = paste0("t = ", t),
      cost = (n - t) * exp(x_grid) - (cumy[n] - cumy[t]) * x_grid + const
    )
  }),
  tibble(
    x = x_grid,
    curve = "global",
    cost = n * exp(x_grid) - cumy[n] * x_grid - beta
  )
)



breaks_curve <- c("global", paste0("t = ", 1:(n - 1)))   # underlying curve names
labels_curve <- paste0("t = ", 0:(n - 1))                # what to display: t = 0,...,t = 19

p2 <- ggplot(cost_df, aes(x = x, y = cost, color = curve)) +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(mini - l, mini + r)) +
  labs(x = "x", y = "Cost(x)", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right") +
  scale_color_discrete(breaks = breaks_curve, labels = labels_curve)

print(p2)
