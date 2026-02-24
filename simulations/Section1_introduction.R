
###
### Section 1. INTRODUCTION
###

#########
# system time, limit = 10s
#########
# BS
# OP
# PELT
# FPOP
# DUST

# install.packages("fpop", repos="http://R-Forge.R-project.org")
library(fpop)
library(fpopw)
library(changepoint)
library(dust)

################################################################################
#####
##### BS system time 10s. How many data points analyzed?
#####


for(i in seq(from = 11.3, to = 11.9, by = 0.1))
{
  n <- i*10^7
  print(n)
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(multiBinSeg(data, Kmax = 9)))
}

###
### response 119.10^6
###

################################################################################
#####
##### OP system time 10s. How many data points analyzed?
#####

for(i in seq(from = 1, to = 1.5, by = 0.1))
{
  n <- i*10^5
  print(n)
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(flat_OP_1D(data,inPenalty = 2*log(n)/2)))
}

###
### response 120 10^3 for 10s max (OK)
###

################################################################################
#####
##### PELT system time 10s. How many data points analyzed?

for(i in seq(from = 2.6, to = 3.2, by = 0.1))
{
  n <- i*10^5
  print(n)
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(cpt.mean(data, method = "PELT", penalty = "BIC")))
}

###
### response 280 10^3 for 10s max (OK)
###

################################################################################
#####
##### FPOP system time 10s. How many data points analyzed?
#####


for(i in seq(from = 3.6, to = 4, by = 0.1))
{
  n <- i*10^7
  print(n)
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(fpop::Fpop(x = data, lambda = 2*log(n))))
}

###
### response 39 10^6
###

################################################################################
#####
##### DUST system time 10s. How many data points analyzed?
#####

for(i in seq(from = 4.2, to = 4.6, by = 0.1))
{
  n <- i*10^7
  print(n)
  cpts <- floor(seq(from = 0.1, to = 1, by = 0.1)*n)
  data <- dataGenerator_1D(chpts = cpts,
                           parameters = c(0,1,0,1,0,1,0,1,0,1),
                           type = "gauss",
                           sdNoise = 1)
  print(system.time(dust.gauss.1D(data, penalty = 2*log(n)/2, model = "gauss")))
}

###
### response 42 10^6
###

# dust.gauss.1D = fast version of dust

