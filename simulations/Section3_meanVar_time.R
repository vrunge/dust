library(dust)
library(parallel)
library(ggplot2)

n <- 10^6
pen <- 2 * log(n)

nb_seg_grid <- c(1, 10, 10^2, 10^3, 10^4)
B <- 100
methods <- c("det_DUST1", "det2_DUST2")

num_cores <- 8

############################

one_run <- function(nb_seg, rep_id, seed)
{
  ## This will now print with PSOCK + outfile=""
  message("nb_seg = ", nb_seg, " | rep = ", rep_id)

  set.seed(seed)
  seg_len <- n / nb_seg
  chpts <- as.integer((1:nb_seg) * seg_len)

  # probabilities of changing at each boundary (tune these)
  p_change <- 0.5
  means <- numeric(nb_seg)
  vars  <- numeric(nb_seg)
  means[1] <- 4
  vars[1]  <- 1

  if(nb_seg > 1)
  {
    for (k in 2:nb_seg)
    {
      means[k] <- means[k-1]
      vars[k]  <- vars[k-1]
      if (runif(1) < p_change){means[k] <- -means[k]}
      else{vars[k]  <- ifelse(vars[k]==1,16,1)}
    }
  }

  y <- dataGenerator_meanVar(
    chpts = chpts,
    means = means,
    sds   = sqrt(vars)
  )

  ord <- sample(methods, length(methods))

  out <- vector("list", length(ord))
  for (k in seq_along(ord))
  {
    m <- ord[k]
    gc(FALSE)

    t <- unname(system.time(res <- dust::dust.meanVar(y, penalty = pen, method = m))[[1]])

    out[[k]] <- data.frame(
      nb_seg = nb_seg,
      rep    = rep_id,
      method = m,
      time_s = t,
      nb_estimated = length(res$changepoints)
    )
  }
  do.call(rbind, out)
}

## ---- Jobs
jobs <- expand.grid(
  nb_seg = nb_seg_grid,
  rep    = seq_len(B),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
base_seed <- 123
jobs$seed <- base_seed + seq_len(nrow(jobs))

#################################################
## Parallel run with PSOCK (prints in RStudio console)

cl <- makeCluster(num_cores, type = "PSOCK", outfile = "")  # <- key for console output

# load needed packages on workers
clusterEvalQ(cl, {
  library(dust)
  NULL
})

# export objects/functions used by workers
clusterExport(
  cl,
  varlist = c("n", "pen", "methods", "jobs", "one_run"),
  envir = environment()
)

res_list <- parLapply(
  cl,
  X = seq_len(nrow(jobs)),
  fun = function(i) one_run(jobs$nb_seg[i], jobs$rep[i], jobs$seed[i])
)

stopCluster(cl)

timings <- do.call(rbind, res_list)


#################################################
## Summaries: time_s + nb_estimated

summ_fun <- function(x)
  c(mean = mean(x), sd = sd(x), median = median(x), min = min(x), max = max(x))

summary_time <- aggregate(
  time_s ~ nb_seg + method,
  data = timings,
  FUN = summ_fun
)

summary_nb <- aggregate(
  nb_estimated ~ nb_seg + method,
  data = timings,
  FUN = summ_fun
)

summary_df_time <- cbind(
  summary_time[, c("nb_seg", "method")],
  as.data.frame(summary_time$time_s)
)

summary_df_nb <- cbind(
  summary_nb[, c("nb_seg", "method")],
  as.data.frame(summary_nb$nb_estimated)
)

summary_df <- merge(
  summary_df_time, summary_df_nb,
  by = c("nb_seg", "method"),
  suffixes = c("_time", "_nb")
)

row.names(summary_df) <- NULL
print(summary_df)

#################################################
## ggplot2 boxplot (with jittered points)

p <- ggplot(timings, aes(x = factor(nb_seg), y = time_s, color = method)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), alpha = 0.2) +
  scale_y_log10() +
  labs(x = "nb_seg", y = "Execution time (s, log10 scale)", color = "Method")

print(p)


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
n <- 10^6
pen <- 2 * log(n)

nb_seg_grid <- c(10^2, 10^3, 10^4)
B <- 100
methods <- c("PELT")

num_cores <- 8

############################
## ---- Jobs
jobs <- expand.grid(
  nb_seg = nb_seg_grid,
  rep    = seq_len(B),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
base_seed <- 123
jobs$seed <- base_seed + seq_len(nrow(jobs))

#################################################
## Parallel run with PSOCK (prints in RStudio console)

cl <- makeCluster(num_cores, type = "PSOCK", outfile = "")  # <- key for console output

# load needed packages on workers
clusterEvalQ(cl, {
  library(dust)
  NULL
})

# export objects/functions used by workers
clusterExport(
  cl,
  varlist = c("n", "pen", "methods", "jobs", "one_run"),
  envir = environment()
)

res_list <- parLapply(
  cl,
  X = seq_len(nrow(jobs)),
  fun = function(i) one_run(jobs$nb_seg[i], jobs$rep[i], jobs$seed[i])
)

stopCluster(cl)

timings <- do.call(rbind, res_list)


#################################################
## Summaries: time_s + nb_estimated

summ_fun <- function(x)
  c(mean = mean(x), sd = sd(x), median = median(x), min = min(x), max = max(x))

summary_time <- aggregate(
  time_s ~ nb_seg + method,
  data = timings,
  FUN = summ_fun
)

summary_nb <- aggregate(
  nb_estimated ~ nb_seg + method,
  data = timings,
  FUN = summ_fun
)

summary_df_time <- cbind(
  summary_time[, c("nb_seg", "method")],
  as.data.frame(summary_time$time_s)
)

summary_df_nb <- cbind(
  summary_nb[, c("nb_seg", "method")],
  as.data.frame(summary_nb$nb_estimated)
)

summary_df <- merge(
  summary_df_time, summary_df_nb,
  by = c("nb_seg", "method"),
  suffixes = c("_time", "_nb")
)

row.names(summary_df) <- NULL
print(summary_df)

#################################################
## ggplot2 boxplot (with jittered points)

p <- ggplot(timings, aes(x = factor(nb_seg), y = time_s, color = method)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), alpha = 0.2) +
  scale_y_log10() +
  labs(x = "nb_seg", y = "Execution time (s, log10 scale)", color = "Method")

print(p)


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################


n <- 10^6
pen <- 4 * log(n)

nb_seg_grid <- c(10)
B <- 64
methods <- c("PELT")

num_cores <- 8

############################
############################
## ---- Jobs
jobs <- expand.grid(
  nb_seg = nb_seg_grid,
  rep    = seq_len(B),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
base_seed <- 123
jobs$seed <- base_seed + seq_len(nrow(jobs))

#################################################
## Parallel run with PSOCK (prints in RStudio console)

cl <- makeCluster(num_cores, type = "PSOCK", outfile = "")  # <- key for console output

# load needed packages on workers
clusterEvalQ(cl, {
  library(dust)
  NULL
})

# export objects/functions used by workers
clusterExport(
  cl,
  varlist = c("n", "pen", "methods", "jobs", "one_run"),
  envir = environment()
)

res_list <- parLapply(
  cl,
  X = seq_len(nrow(jobs)),
  fun = function(i) one_run(jobs$nb_seg[i], jobs$rep[i], jobs$seed[i])
)

stopCluster(cl)

timings <- do.call(rbind, res_list)


#################################################
## Summaries: time_s + nb_estimated

summ_fun <- function(x)
  c(mean = mean(x), sd = sd(x), median = median(x), min = min(x), max = max(x))

summary_time <- aggregate(
  time_s ~ nb_seg + method,
  data = timings,
  FUN = summ_fun
)

summary_nb <- aggregate(
  nb_estimated ~ nb_seg + method,
  data = timings,
  FUN = summ_fun
)

summary_df_time <- cbind(
  summary_time[, c("nb_seg", "method")],
  as.data.frame(summary_time$time_s)
)

summary_df_nb <- cbind(
  summary_nb[, c("nb_seg", "method")],
  as.data.frame(summary_nb$nb_estimated)
)

summary_df <- merge(
  summary_df_time, summary_df_nb,
  by = c("nb_seg", "method"),
  suffixes = c("_time", "_nb")
)

row.names(summary_df) <- NULL
print(summary_df)

#################################################
## ggplot2 boxplot (with jittered points)

p <- ggplot(timings, aes(x = factor(nb_seg), y = time_s, color = method)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.6) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), alpha = 0.2) +
  scale_y_log10() +
  labs(x = "nb_seg", y = "Execution time (s, log10 scale)", color = "Method")

print(p)

