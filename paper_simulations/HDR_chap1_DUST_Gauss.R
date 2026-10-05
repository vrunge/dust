
###
### HDR Chapter 1 -- DUST, Gaussian model
###
### Back to dust0's original design (see
### /Users/vrunge/Library/CloudStorage/Dropbox/A_HDR/algo_DUST/dust0):
### a manual linear scan, not the bisection search / curve-fit machinery
### built up in Section1_introduction_gauss.R. For each (algorithm,
### nb_seg) cell below: run the for loop, read off which n the printed
### system.time() crosses 10s at, write it into the "### response" line.
### One cell at a time, by hand -- this file is NOT meant to be run
### start to finish unattended.
###
### Same 3 experiments as Section1_introduction_gauss.R (10, 100, 1000
### equal segments, means alternating 0/1, sdNoise=1, penalty=2*log(n)).
### OP and PELT each get two blocks -- dust's own engine, and the
### changepoint package -- so both can be compared the same way the
### Poisson file compares them.
###
### Ranges below are seeded from this session's own bisection/curve-fit
### results where available (dust BS/OP/PELT/FPOP/DUST -- see
### Section1_introduction_gauss.R's trailing comment block for the exact
### figures) and narrowed to a human-scannable handful of steps around
### them. The changepoint-package OP/PELT ranges were NEVER measured for
### Gaussian in this session (only for Poisson) -- they are rough,
### labelled guesses to start the manual search from, nothing more.

library(fpop)
library(changepoint)
library(dust)

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  means <- rep(c(0, 1), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = means, type = "gauss", sdNoise = 1)
}

################################################################################
################################################################################
##### nb_seg = 10
################################################################################
################################################################################

nb_seg <- 10

################################################################################
##### BS (fpop::multiBinSeg) -- anchor: 172015625 (measured)
################################################################################


for (i in seq(from = 17, to = 19, by = 0.5))
{
  n <- i * 10^7
  print(n)
  temp <- 0
  for(j in 1:10)
  {
    data <- build_data(n, nb_seg)
    temp <- temp + system.time(fpop::multiBinSeg(data, Kmax = nb_seg))
    print(j)
  }
  print(c(i,temp/10))
}

### response: ___


################################################################################
##### OP, dust (dust::dust.1D, method="OP", backend="scalar") -- anchor: 285957 (measured)
################################################################################

for (i in seq(from = 2.6, to = 3.1, by = 0.1))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, penalty = 2*log(n), model = "gauss", method = "OP", backend = "scalar")))
}

### response: ___


################################################################################
##### OP, changepoint (cpt.mean, method="SegNeigh") -- UNCALIBRATED GUESS
################################################################################

for (i in seq(from = 1.0, to = 2.5, by = 0.25))
{
  n <- i * 10^4
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.mean(data, method = "SegNeigh", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### PELT, dust (dust::dust.1D, method="PELT", backend="scalar") -- anchor: 285957 (measured)
################################################################################

for (i in seq(from = 2.6, to = 3.1, by = 0.1))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, penalty = 2*log(n), model = "gauss", method = "PELT", backend = "scalar")))
}

### response: ___


################################################################################
##### PELT, changepoint (cpt.mean, method="PELT") -- UNCALIBRATED GUESS
################################################################################

for (i in seq(from = 2.0, to = 3.5, by = 0.25))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.mean(data, method = "PELT", penalty = "BIC")))
}

### response: ___


################################################################################
##### FPOP (fpop::Fpop) -- anchor: 56687500 (measured)
################################################################################

for (i in seq(from = 5.4, to = 5.95, by = 0.1))
{
  n <- i * 10^7
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(fpop::Fpop(x = data, lambda = 2*log(n))))
}

### response: ___


################################################################################
##### DUST (dust::dust.1D, model="gauss", method="DUST", backend = "highway") -- anchor: 83648438 (measured)
################################################################################

for (i in seq(from = 8.0, to = 8.7, by = 0.1))
{
  n <- i * 10^7
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, 2*log(n), model = "gauss", method = "DUST", backend = "highway")))
}

### response: ___



################################################################################
################################################################################
##### nb_seg = 100
################################################################################
################################################################################

nb_seg <- 100

################################################################################
##### BS (fpop::multiBinSeg) -- anchor: 37146484 (measured)
################################################################################

for (i in seq(from = 3.5, to = 3.95, by = 0.05))
{
  n <- i * 10^7
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(fpop::multiBinSeg(data, Kmax = nb_seg)))
}

### response: ___


################################################################################
##### OP, dust -- anchor: 911719 (measured)
################################################################################

for (i in seq(from = 8.6, to = 9.6, by = 0.1))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, penalty = 2*log(n), model = "gauss", method = "OP", backend = "scalar")))
}

### response: ___


################################################################################
##### OP, changepoint (SegNeigh) -- UNCALIBRATED GUESS
################################################################################

for (i in seq(from = 3.0, to = 7.0, by = 1.0))
{
  n <- i * 10^4
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.mean(data, method = "SegNeigh", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### PELT, dust -- anchor: 911719 (measured)
################################################################################

for (i in seq(from = 8.6, to = 9.6, by = 0.1))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, penalty = 2*log(n), model = "gauss", method = "PELT", backend = "scalar")))
}

### response: ___


################################################################################
##### PELT, changepoint -- UNCALIBRATED GUESS
################################################################################

for (i in seq(from = 7.0, to = 11.0, by = 0.5))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.mean(data, method = "PELT", penalty = "BIC")))
}

### response: ___


################################################################################
##### FPOP -- anchor: 62101563 (measured)
################################################################################

for (i in seq(from = 5.9, to = 6.5, by = 0.1))
{
  n <- i * 10^7
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(fpop::Fpop(x = data, lambda = 2*log(n))))
}

### response: ___


################################################################################
##### DUST -- anchor: 92960938 (measured)
################################################################################

for (i in seq(from = 8.9, to = 9.6, by = 0.1))
{
  n <- i * 10^7
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, 2*log(n), model = "gauss", method = "DUST", backend = "highway")))
}

### response: ___



################################################################################
################################################################################
##### nb_seg = 1000
################################################################################
################################################################################

nb_seg <- 1000

################################################################################
##### BS (fpop::multiBinSeg) -- anchor: 4012207 (measured)
################################################################################

for (i in seq(from = 3.7, to = 4.3, by = 0.1))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(fpop::multiBinSeg(data, Kmax = nb_seg)))
}

### response: ___


################################################################################
##### OP, dust -- anchor: 2343750 (measured)
################################################################################

for (i in seq(from = 2.2, to = 2.5, by = 0.05))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, penalty = 2*log(n), model = "gauss", method = "OP", backend = "scalar")))
}

### response: ___


################################################################################
##### OP, changepoint (SegNeigh) -- UNCALIBRATED GUESS
################################################################################

for (i in seq(from = 1.0, to = 2.0, by = 0.2))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.mean(data, method = "SegNeigh", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### PELT, dust -- anchor: 2750000 (measured)
################################################################################

for (i in seq(from = 2.6, to = 2.9, by = 0.05))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, penalty = 2*log(n), model = "gauss", method = "PELT", backend = "scalar")))
}

### response: ___


################################################################################
##### PELT, changepoint -- UNCALIBRATED GUESS
################################################################################

for (i in seq(from = 1.8, to = 3.0, by = 0.2))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.mean(data, method = "PELT", penalty = "BIC")))
}

### response: ___


################################################################################
##### FPOP -- anchor: 72156250 (measured)
################################################################################

for (i in seq(from = 6.9, to = 7.5, by = 0.1))
{
  n <- i * 10^7
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(fpop::Fpop(x = data, lambda = 2*log(n))))
}

### response: ___


################################################################################
##### DUST -- anchor: 106929688 (measured)
################################################################################

for (i in seq(from = 10.2, to = 11.2, by = 0.1))
{
  n <- i * 10^7
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(dust.1D(data, 2*log(n), model = "gauss", method = "DUST", backend = "highway")))
}

### response: ___
