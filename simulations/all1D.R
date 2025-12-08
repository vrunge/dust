





#####################################################################

devtools::install_github("vrunge/dust")
library(dust)
n <- 10^4
pen <- 2*log(n)
data <- dataGenerator_1D(chpts = n,
                         parameters = 1,
                         type = "gauss")
system.time(dust.1D(data, penalty = pen, method = "OP"))
# 0.5 s


######################################################################


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
                              means = sample(c(0,0), 50, replace = TRUE),
                              sds = sample(c(1,1), 50, replace = TRUE))
plot(data, type = 'l')
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




