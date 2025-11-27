

for(mod in c("gauss","poisson","exp","geom","bern","binom","negbin","variance"))
{
  print(c())
  print(c())
  print(rep(mod,15))
  print(n)

  n <- 10^4
  pen <- 2*log(n)/2
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

  res0r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTr")
  res1r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUST")
  res2r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTgs")
  res3r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTqn")
  res4r <- dust.1D(data, penalty = pen, model = mod, method = "rand_PELT")
  res5r <- dust.1D(data, penalty = pen, model = mod, method = "rand_DUSTbs")

  print(c(all(res0d$changepoints == res1d$changepoints),
  all(res0d$changepoints == res2d$changepoints),
  all(res0d$changepoints == res3d$changepoints),
  all(res0d$changepoints == res4d$changepoints),
  all(res0d$changepoints == res5d$changepoints)))

  print(c(all(res0r$changepoints == res1r$changepoints),
  all(res0r$changepoints == res2r$changepoints),
  all(res0r$changepoints == res3r$changepoints),
  all(res0r$changepoints == res4r$changepoints),
  all(res0r$changepoints == res5r$changepoints)))
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




