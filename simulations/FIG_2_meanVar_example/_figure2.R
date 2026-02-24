
library(dust)

n <- 10^2
beta <- 2*log(n)
y <-  dataGenerator_meanVar(chpts = c(n))
y <- c(rep(0, n), rep(10, n))
plot(y)

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

### DUST : randIndex_Eval0
### PELT : randIndex_Eval5


############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

oneSimu <- function(n)
{
  beta <- 4*log(n)
  y <-  dataGenerator_meanVar(chpts = n)
  res <- dust.meanVar(data = y, penalty = beta, method = "det_DUST1")
  cat("0")
  return(res$nb)
}


oneSimu2 <- function(n)
{
  beta <- 4*log(n)
  y <-  dataGenerator_meanVar(chpts = c(round(1/3*n), round(2/3*n),n),
                              means = c(0, 0, 0.5),
                              sds = c(1, 2, 2))
  res <- dust.meanVar(data = y, penalty = beta, method = "det_DUST1")
  cat("0" )
  return(res$nb)
}


# Load required libraries
library(ggplot2)

# Parameters
n_steps <- 10^4  # Number of time steps
n_simulations <- 100  # Number of random walks


#####################################################################

dust_simus <- replicate(n_simulations, oneSimu(n_steps))

beta <- 4*log(n_steps)
y <-  dataGenerator_meanVar(chpts = n_steps)
res <-  dust::dust.meanVar(data = y, penalty = beta, method = "DUST1")


# Calculate the mean and confidence intervals for each time step
means <- apply(dust_simus, 1, mean)
up <- apply(dust_simus, 1, function(x) quantile(x,0.95))
down <- apply(dust_simus, 1, function(x) quantile(x,0.05))

# Create a data frame for plotting
time <- 1:n_steps
df <- data.frame(time = time,
                 examplarOP = c(res$nb),  # Pick the first simulation as the exemplar
                 examplar = dust_simus[,1],  # Pick the first simulation as the exemplar
                 mean = means,
                 lower = down,
                 upper = up)
y_min <- min(c(df$lower, df$examplar))  # You can set a custom value if needed, e.g., -20
y_max <- max(c(df$upper, df$examplar)) # You can set a custom value if needed, e.g., 20


ggplot() +
  geom_line(data = df, aes(x = time, y = examplarOP, color = "PELT"), size = 1) +  # Exemplar random walk
  geom_line(data = df, aes(x = time, y = examplar, color = "one dust"), size = 1) +  # Exemplar random walk
  geom_line(data = df, aes(x = time, y = mean, color = "1000 dust mean"), size = 1) +  # Mean curve
  geom_ribbon(data = df,aes(x = time, ymin = lower, ymax = upper), fill = "grey60", alpha = 0.4) +  # Confidence interval
  labs(title = "",
       x = "",
       y = "number of indices",
       caption = "",
       color = "") + ylim(y_min, y_max) +
  scale_color_manual(values = c("PELT" = "orange", "one dust" = "blue", "1000 dust mean" = "red"),  breaks = c("PELT", "one dust", "1000 dust mean")) + # Custom colors for the legend
  theme(
    axis.text = element_text(size = 16),            # Increase size of axis numbers
    axis.title.y = element_text(size = 20),         # Increase size of the x-axis label
    legend.text = element_text(size = 20),          # Increase size of the legend text
    legend.title = element_text(size = 16),          # Increase size of the legend title
    legend.position = "top"
    )

##################
##########################
###########################################
##########################
##################


dust_simusB <- replicate(n_simulations, oneSimu2(n_steps))

beta <- 4*log(n_steps)
yB <-  dataGenerator_meanVar(chpts = c(round(1/3*n_steps), round(2/3*n_steps),n_steps),
                            means = c(0, 0, 0.5),
                            sds = c(1, 2, 2))
resB <- dust.partitioner.meanVar(method = "randIndex_Eval5")$quick(data = yB, penalty = beta)


# Calculate the mean and confidence intervals for each time step
meansB <- apply(dust_simusB, 1, mean)
upB <- apply(dust_simusB, 1, function(x) quantile(x,0.95))
downB <- apply(dust_simusB, 1, function(x) quantile(x,0.05))


# Create a data frame for plotting
time <- 1:n_steps
dfB <- data.frame(time = time,
                 examplarOP = c(resB$nb),  # Pick the first simulation as the exemplar
                 examplar = dust_simusB[,1],  # Pick the first simulation as the exemplar
                 mean = meansB,
                 lower = downB,
                 upper = upB)
y_min <- min(dfB$lower)  # You can set a custom value if needed, e.g., -20
y_max <- max(dfB$upper)  # You can set a custom value if needed, e.g., 20


ggplot() +
  geom_line(data = dfB, aes(x = time, y = examplarOP, color = "PELT"), size = 1) +  # Exemplar random walk # Exemplar random walk
  geom_line(data = dfB, aes(x = time, y = mean, color = "mean"), size = 1) +  # Mean curve
  geom_ribbon(data = dfB, aes(x = time, ymin = lower, ymax = upper), fill = "grey60", alpha = 0.5) +  # Confidence interval
  labs(title = "",
       x = "",
       y = "number of indices",
       caption = "",
       color = "",
       size = 20) +
  scale_color_manual(values = c("PELT" = "orange", "mean" = "red"),  breaks = c("PELT", "one dust", "mean"))+  # Custom colors for the legend
  theme(
    axis.text = element_text(size = 16),            # Increase size of axis numbers
    axis.title.y = element_text(size = 20),         # Increase size of the x-axis label
    legend.text = element_text(size = 20),          # Increase size of the legend text
    legend.title = element_text(size = 16),          # Increase size of the legend title
    legend.position = "top"
  ) +
  geom_vline(xintercept = resB$changepoints[1], linetype = "twodash", size = 0.5) +  # Add a vertical line at x = 2
  geom_vline(xintercept = resB$changepoints[2], linetype = "twodash", size = 0.5) +  # Add another vertical line at x = 3
  geom_vline(xintercept = resB$changepoints[3], linetype = "twodash", size = 0.5)   # Add anot


resB$changepoints
plot(yB)
plot(resB$nb)


means[n_steps]/n_steps * 100
n_steps*(n_steps+1)/2/sum(means)
n_steps*(n_steps+1)/2/sum(meansB)




