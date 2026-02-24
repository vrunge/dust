library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------
# Data generators
# -----------------------------

data_gauss_noPruning <- function(n, beta)
{
  t <- seq_len(n)
  sqrt(beta / n) * (sqrt(n - 1) - sqrt(t * (n - t)) + sqrt((t - 1) * (n - t + 1)))
}

data_exponential_noPruning <- function(n, beta, Y)
{
  # For each t=1,...,n-1 solve in x in (0, Y):
  # -log((Y - (t/n)*x)/(Y - (t/n)*Y)) + (t/n)*log((Y - (t/n)*x)/(x - (t/n)*x)) - beta/n = 0

  f_factory <- function(t)
  {
    function(x)
    {
      tn <- t / n
      # Terms rewritten for clarity; keep original math
      -log((Y - tn * x) / (Y - tn * Y)) +
        tn * log((Y - tn * x) / (x - tn * x)) -
        beta / n
    }
  }

  y <- numeric(n)
  for (t in 1:(n - 1))
  {
    # uniroot requires a sign change; keep same interval as your code
    y[t] <- uniroot(f_factory(t), interval = c(0, Y))$root
  }
  y[n] <- Y

  # Your original post-processing:
  # res = diff((1:n)*y), then prepend y[1]
  res <- c(y[1], diff(seq_len(n) * y))
  res
}

# -----------------------------
# Parameters + compute series
# -----------------------------

n <- 20
Y <- 1
beta_gauss <- n/(n-1) * Y^2

d1 <- data_gauss_noPruning(n, beta = beta_gauss)
d2 <- data_exponential_noPruning(n, beta = 0.1, Y = Y)
d3 <- data_exponential_noPruning(n, beta = 1,   Y = Y)

# Means (as in your script)
mean(d1); mean(d2); mean(d3)

# -----------------------------
# Tidy data for ggplot2
# -----------------------------

df <- tibble(
  t = seq_len(n),
  Gauss = d1,
  `Exponential, beta = 0.1` = d2,
  `Exponential, beta = 1`   = d3
) |>
  pivot_longer(-t, names_to = "model", values_to = "value")

# -----------------------------
# Plot (lines + points like type = "b")
# -----------------------------

ggplot(df, aes(x = t, y = value, color = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.6) +
  scale_x_continuous(breaks = seq_len(n)) +
  labs(
    x = "t",
    y = "value",
    color = NULL
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###### ### ### ### ### ### ### ###

### verification

### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###


pen <- 2
resG <- data_gauss_noPruning(10, pen)
resE <- data_exponential_noPruning(10, pen, 1)

plot(resG)

### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ###

### REMARK
### we consider "-theta" instead of "theta"
### the curves are in decreasing order
###

plot_data_exponential_noPruning <- function(res)
{
  n <- length(res)
  cumy <- cumsum(res)
  u <- cumy[1]
  res <- 1 + log(u)
  x <- 1/((cumy[n] - cumy[1])/(n-1))
  mini <- (n-1)*(-log(x)) + (cumy[n] - cumy[1])*x + res


  for(t in 1:(n-1))
  {
    Cost <- function(x)
    {
      u <- cumy[t]/t
      res <- t*(1 + 1 *log(u))
      return(-(n-t)*log(x) + (cumy[n] - cumy[t])*x + res)
    }

    curve(Cost, n = 1000, col = t+1, from = 0, to = 1, ylim = c(mini-0.01, mini+3))
    par(new = TRUE)
  }
  Cost <- function(x)
  {
    return(-(n)*log(x) + (cumy[n])*x - pen)
  }
  curve(Cost, n = 1000, col = 1, from = 0, to = 1, ylim = c(mini-0.01, mini+3))


  par(new = FALSE)
}


plot_data_exponential_noPruning(resE)

