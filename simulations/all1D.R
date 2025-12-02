obj1 <- new(DUST_meanVar, "det_DUST1")
obj1$get_info()

obj2 <- new(DUST_meanVar, "det2_DUST2")
obj2$get_info()

### MEAN AND VAR
mod = "gauss"
n <- 10^6
pen <- 2*log(n)
cpts <- n
data <- dataGenerator_meanVar(chpts = 1:50*(n/50),
                              means = sample(c(0,1), 50, replace = TRUE),
                              sds = sample(c(1,1.5), 50, replace = TRUE))

#res1d <- dust.meanVar(data, penalty = pen, method = "PELT")
#res2d <- dust.meanVar(data, penalty = pen, method = "det_DUSTr") #red
#res3d <- dust.meanVar(data, penalty = pen, method = "det_DUST") #green
system.time(res4d <- dust.meanVar(data, penalty = pen, method = "det_DUST1")) #blue
system.time(res5d <- dust.meanVar(data, penalty = pen, method = "det2_DUST2")) #blue

res4d$changepoints

#plot(res1d$nb, col = 1, type = 'l', ylim = c(0, max(res2d$nb, res3d$nb)))
#par(new = TRUE)
MAX <- max(res4d$nb, res5d$nb)

#plot(res2d$nb, col = 2, type = 'l', ylim = c(0, MAX))
#par(new = TRUE)

plot(res4d$nb, col = 3, type = 'l', ylim = c(0, MAX), main = "number of non-pruned indices over time")
par(new = TRUE)
plot(res5d$nb, col = 4, type = 'l', ylim = c(0, MAX))
legend("topleft",
       legend = c("1D dual (0.05% indices left. 10s) ", "2D dual  (0.02% indices left. 10s)"),  # change labels as you like
       col    = c(3, 4),
       lty    = 1,
       lwd = 2,
       bty    = "n")

res4d$nb[n]/n
res5d$nb[n]/n

length(res4d$changepoints)

#sum(abs(res3d$costQ - res2d$costQ)[-1])
#sum(abs(res1d$costQ - res2d$costQ)[-1])

res4d$changepoints
res5d$changepoints
all(res4d$changepoints == res5d$changepoints)
res4d$nb[n]/n
res5d$nb[n]/n


res5d$costQ


all(res2d$changepoints == res3d$changepoints)
all(res3d$changepoints == res4d$changepoints)



res4d$nb - res2d$nb


#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################



res4d


plot(data, type = 'l')

plot(res0d$nb, type = 'l', ylim = c(0,max(res0d$nb)), main = mod)
plot(res1d$nb, type = 'l', ylim = c(0,max(res1d$nb)), main = mod)







for(mod in c("gauss","poisson","exp","geom","bern","binom","negbin","variance"))
{
  print(c())
  print(c())
  print(rep(mod,15))
  n <-10^4
  pen <- 2*log(n)
  cpts <- n
  #cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = 0.5,
                           type = mod, nbTrials = 5)
  data <- data_normalization_1D(data, type = mod)
  #print(data)
  #res0d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTr")
  #res1d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUST")
  res2d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTgs")
  res3d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTqn")
  #res4d <- dust.1D(data, penalty = pen, model = mod, method = "det_PELT")
  res5d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTbs")
  res6d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTib")

  print(c(#all(res0d$changepoints == res1d$changepoints),
          all(res3d$changepoints == res2d$changepoints),
          #all(res0d$changepoints == res3d$changepoints),
          #all(res0d$changepoints == res4d$changepoints),
          all(res3d$changepoints == res5d$changepoints),
          all(res3d$changepoints == res6d$changepoints)))
  print(res3d$changepoints)

  print(sum(abs(res3d$nb - res6d$nb)))
  print(sum(res3d$nb - res6d$nb))

  plot(res3d$nb, type = 'l', ylim = c(0,max(res3d$nb,res6d$nb)), main = mod)
  par(new = TRUE)
  plot(res6d$nb, col = 2, type = 'l', pch = '-', ylim = c(0,max(res3d$nb,res6d$nb)))
}

res3d$costQ - res6d$costQ
sum(abs(res3d$nb - res6d$nb))
res3d$changepoints
res6d$changepoints



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


all_times <- list()  # optional: to store timings per model

for (mod in c("gauss","poisson","exp","geom","bern","binom","negbin","variance"))
#for (mod in c("variance"))
{
  cat("\n============================\n")
  cat("Model:", mod, "\n")

  n   <- 1*10^6
  pen <- 2 * log(n)
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1) * n)

  data <- dataGenerator_1D(
    chpts      = cpts,
    parameters = 0.5 * c(1,1,1,1,1,1,1,1,1,1),
    type       = mod
  )
  data <- data_normalization_1D(data, type = mod)

  ##------------------##
  ## deterministic    ##
  ##------------------##
  #t0d <- system.time(res0d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTr"))
  #t1d <- system.time(res1d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUST"))
  #t2d <- system.time(res2d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTgs"))
  t3d <- system.time(res3d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTqn"))
  t5d <- system.time(res5d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTbs"))
  t6d <- system.time(res6d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTib"))

  times_det <- rbind(
    #det_DUSTr  = t0d,
    #det_DUST   = t1d,
    #det_DUSTgs = t2d,
    det_DUSTqn = t3d,
    det_DUSTbs = t5d,
    det_DUSTib = t6d
  )

  ##------------------##
  ## randomized       ##
  ##------------------##
  #t0r <- system.time(res0r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTr"))
  #t1r <- system.time(res1r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUST"))
  #t2r <- system.time(res2r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTgs"))
  #t3r <- system.time(res3r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTqn"))
  #t5r <- system.time(res5r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTbs"))
  #t6r <- system.time(res6r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTib"))

  times_rand <- rbind(
    #rand_DUSTr  = t0r,
   # rand_DUST   = t1r,
    #rand_DUSTgs = t2r,
   # rand_DUSTqn = t3r,
    #rand_DUSTbs = t5r,
    #rand_DUSTib = t6r
  )

  ##------------------##
  ## print results    ##
  ##------------------##
  cat("\nDeterministic timings (user/sys/elapsed):\n")
  print(times_det)

  #cat("\nRandomized timings (user/sys/elapsed):\n")
  #print(times_rand)

  cat("\nDeterministic changepoint agreement (vs det_DUSTr):\n")
  print(c(
    #all(res0d$changepoints == res1d$changepoints),
    #all(res0d$changepoints == res2d$changepoints),
    #all(res0d$changepoints == res3d$changepoints),
    #all(res0d$changepoints == res5d$changepoints)
    #all(res0d$changepoints == res6d$changepoints)
  ))

  #cat("\nRandomized changepoint agreement (vs rand_DUSTr):\n")
  #print(c(
    #all(res0r$changepoints == res1r$changepoints),
    #all(res0r$changepoints == res2r$changepoints),
    #all(res0r$changepoints == res3r$changepoints),
    #all(res0r$changepoints == res5r$changepoints)
    #all(res0r$changepoints == res6r$changepoints)
  #))

  #cat("\nRandomized changepoints (rand_DUSTr):\n")
  print(res3d$changepoints)

  ## optional: store timings in a list for later analysis
  #all_times[[mod]] <- list(det = times_det, rand = times_rand)
  all_times[[mod]] <- list(det = times_det)
}


all_times





################################################################################
################################################################################
################################################################################
################################################################################




for(mod in c("gauss","poisson","exp","geom","bern","binom","negbin","variance"))
{
  print(c())
  print(c())
  print(rep(mod,15))
  print(n)

  n <- 10^4
  pen <- 2*log(n)/3
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = 0.5*c(0.1,1,0.1,1,0.1,1,0.1,1,0.1,1),
                           type = mod)
  data <- data_normalization_1D(data, type = mod)
  res0d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTr")
  res1d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUST")
  res2d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTgs")
  res3d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTqn")
  res4d <- dust.1D(data, penalty = pen, model = mod, method = "det_PELT")
  res5d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTbs")
  res6d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTib")

  res0r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTr")
  res1r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUST")
  res2r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTgs")
  res3r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTqn")
  res4r <- dust.1D(data, penalty = pen, model = mod, method = "rand_PELT")
  res5r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTbs")
  res6r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTib")

  print(c(all(res0d$changepoints == res1d$changepoints),
  all(res0d$changepoints == res2d$changepoints),
  all(res0d$changepoints == res3d$changepoints),
  all(res0d$changepoints == res4d$changepoints),
  all(res0d$changepoints == res5d$changepoints),
  all(res0d$changepoints == res6d$changepoints)))

  print(c(all(res0r$changepoints == res1r$changepoints),
  all(res0r$changepoints == res2r$changepoints),
  all(res0r$changepoints == res3r$changepoints),
  all(res0r$changepoints == res4r$changepoints),
  all(res0r$changepoints == res5r$changepoints),
  all(res0r$changepoints == res6r$changepoints)))
  print(res0r$changepoints)
}














obj_dust <- dust.object.1D(
  model = "gfaussfn",
  method = "det_DUSTgs",
  nbLoops = 0)
obj_dust$get_info()



obj_dust <- dust.object.1D(
  model = "gauss",
  method = "PELT",
  nbLoops = 0)
obj_dust$get_info()


obj_dust <- dust.object.1D(
  model = "gauss",
  method = "OP",
  nbLoops = 0)
obj_dust$get_info()







mod  <- c("gauss")
n <- 10^3
  print(n)
  pen <- 2*log(n)
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = 5*c(0.1,1,0.1,1,0.1,1,0.1,1,0.1,1),
                           type = mod)
  data <- data_normalization_1D(data, type = mod)
  res0d <- dust.1D(data, penalty = pen, model = mod, method = "det_DUSTr")


  plot(data, type = 'l')
  abline(v = res0d$changepoints, col = "red")




