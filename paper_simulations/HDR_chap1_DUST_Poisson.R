
###
### HDR Chapter 1 -- DUST, Poisson model
###
### Back to dust0's original design (see
### /Users/vrunge/Library/CloudStorage/Dropbox/A_HDR/algo_DUST/dust0):
### a manual linear scan, not the bisection search / curve-fit machinery
### built up in Section1_introduction_poisson.R. For each (algorithm,
### nb_seg) cell below: run the for loop, read off which n the printed
### system.time() crosses 10s at, write it into the "### response" line.
### One cell at a time, by hand -- this file is NOT meant to be run
### start to finish unattended.
###
### Same 3 experiments as Section1_introduction_poisson.R (10, 100, 1000
### equal segments, rates alternating 2/8, penalty=2*log(n)). OP and
### PELT each get two blocks -- dust's own engine, and the changepoint
### package. GFPOP (gfpop package) stands in for FPOP, which is
### Gaussian-only. DUST uses data_normalization_1D() first, matching
### dust's own convention; BS/OP/PELT/GFPOP get the raw Poisson counts.
###
### Ranges below:
###  - nb_seg=10: seeded from this session's own measured curve-fit
###    results (see Section1_introduction_poisson.R / conversation),
###    narrowed to a human-scannable handful of steps. Also verified
###    this session that every method here detects the true changepoints
###    correctly at these n (not just "finishes fast" -- see conversation).
###  - nb_seg=100 and nb_seg=1000: NEVER MEASURED for Poisson this
###    session. Ranges are extrapolated from the Gaussian experiment's
###    OWN measured growth pattern between segment counts (BS shrinks
###    sharply; OP/PELT/GFPOP/DUST grow moderately) applied to the
###    nb_seg=10 Poisson anchors. Treat these as a starting point for
###    the manual search, not a prediction to trust -- Poisson's
###    growth pattern was never confirmed to match Gaussian's.
###  - BS is crash-ceiling-limited at nb_seg=10 (changepoint's BinSeg
###    hits a C stack overflow well before 10s of actual compute --
###    see conversation), not time-limited: the printed system.time()
###    values will likely be far under 10s right up to where it crashes.
###    Expect an R error, not a slow run, if the range is too high.

library(gfpop)
library(changepoint)
library(dust)

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  rates <- rep(c(4, 6), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = rates, type = "poisson")
}

################################################################################
################################################################################
##### nb_seg = 10
################################################################################
################################################################################

nb_seg <- 10

################################################################################
##### BS (changepoint::cpt.meanvar, method="BinSeg") -- anchor: 984609
##### (CRASH ceiling, not time ceiling -- see note above)
################################################################################

for (i in seq(from = 9.0, to = 10.0, by = 0.2))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "BinSeg", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### OP, dust (dust.1D, model="poisson", method="OP", backend = "highway") -- anchor: 172525 (measured)
################################################################################

for (i in seq(from = 1.5, to = 2.0, by = 0.05))
{
  n <- i * 10^5
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "OP", backend = "highway")))
}

### response: ___


################################################################################
##### OP, changepoint (SegNeigh) -- anchor: 7805 (measured)
################################################################################

for (i in seq(from = 6.0, to = 9.0, by = 0.5))
{
  n <- i * 10^3
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "SegNeigh", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### PELT, dust -- anchor: 190767 (measured)
################################################################################

for (i in seq(from = 1.7, to = 2.1, by = 0.05))
{
  n <- i * 10^5
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "PELT", backend = "highway")))
}

### response: ___


################################################################################
##### PELT, changepoint -- anchor: 203897 (measured)
################################################################################

for (i in seq(from = 1.8, to = 2.3, by = 0.05))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "PELT", penalty = "BIC")))
}

### response: ___


################################################################################
##### GFPOP (gfpop::gfpop, type="poisson") -- anchor: 2782385 (measured)
################################################################################

for (i in seq(from = 2.6, to = 3.0, by = 0.05))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  g <- graph(type = "std", penalty = 2*log(n))
  print(system.time(gfpop(data, mygraph = g, type = "poisson")))
}

### response: ___


################################################################################
##### DUST -- anchor: 24168510 (measured)
################################################################################

for (i in seq(from = 2.3, to = 2.55, by = 0.05))
{
  n <- i * 10^7
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "DUST", backend = "highway")))
}

### response: ___



################################################################################
################################################################################
##### nb_seg = 100  -- ranges EXTRAPOLATED from the Gaussian 10->100 growth
##### pattern applied to the nb_seg=10 Poisson anchors above. UNMEASURED.
################################################################################
################################################################################

nb_seg <- 100

################################################################################
##### BS -- crash-ceiling behaviour; unknown if/how it scales with nb_seg.
##### Range is a broad guess, not a narrowed one like the others.
################################################################################

for (i in seq(from = 2.0, to = 10.0, by = 1.0))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "BinSeg", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### OP, dust -- extrapolated guess ~550000 (x3.19 vs nb_seg=10, Gaussian factor)
################################################################################

for (i in seq(from = 4.5, to = 6.5, by = 0.25))
{
  n <- i * 10^5
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "OP", backend = "highway")))
}

### response: ___


################################################################################
##### OP, changepoint (SegNeigh) -- extrapolated guess ~25000
################################################################################

for (i in seq(from = 1.5, to = 3.5, by = 0.25))
{
  n <- i * 10^4
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "SegNeigh", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### PELT, dust -- extrapolated guess ~608000
################################################################################

for (i in seq(from = 5.0, to = 7.0, by = 0.25))
{
  n <- i * 10^5
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "PELT", backend = "highway")))
}

### response: ___


################################################################################
##### PELT, changepoint -- extrapolated guess ~650000
################################################################################

for (i in seq(from = 5.0, to = 7.5, by = 0.25))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "PELT", penalty = "BIC")))
}

### response: ___


################################################################################
##### GFPOP -- extrapolated guess ~3050000
################################################################################

for (i in seq(from = 2.8, to = 3.3, by = 0.05))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  g <- graph(type = "std", penalty = 2*log(n))
  print(system.time(gfpop(data, mygraph = g, type = "poisson")))
}

### response: ___


################################################################################
##### DUST -- extrapolated guess ~26860000
################################################################################

for (i in seq(from = 2.5, to = 2.9, by = 0.05))
{
  n <- i * 10^7
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "DUST", backend = "highway")))
}

### response: ___



################################################################################
################################################################################
##### nb_seg = 1000  -- ranges EXTRAPOLATED further (100->1000 Gaussian
##### factor applied on top of the nb_seg=100 guesses above). UNMEASURED.
################################################################################
################################################################################

nb_seg <- 1000

################################################################################
##### BS -- crash-ceiling behaviour; unknown if/how it scales with nb_seg.
##### Range is a broad guess.
################################################################################

for (i in seq(from = 0.5, to = 5.0, by = 0.5))
{
  n <- i * 10^5
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "BinSeg", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### OP, dust -- extrapolated guess ~1414000
################################################################################

for (i in seq(from = 1.2, to = 1.7, by = 0.05))
{
  n <- i * 10^6
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "OP", backend = "highway")))
}

### response: ___


################################################################################
##### OP, changepoint (SegNeigh) -- extrapolated guess ~64000
################################################################################

for (i in seq(from = 4.0, to = 9.0, by = 0.5))
{
  n <- i * 10^4
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "SegNeigh", penalty = "BIC", Q = 2*nb_seg + 10)))
}

### response: ___


################################################################################
##### PELT, dust -- extrapolated guess ~1837000
################################################################################

for (i in seq(from = 1.6, to = 2.1, by = 0.05))
{
  n <- i * 10^6
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "PELT", backend = "highway")))
}

### response: ___


################################################################################
##### PELT, changepoint -- extrapolated guess ~1964000
################################################################################

for (i in seq(from = 1.7, to = 2.3, by = 0.05))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  print(system.time(cpt.meanvar(data, test.stat = "Poisson", method = "PELT", penalty = "BIC")))
}

### response: ___


################################################################################
##### GFPOP -- extrapolated guess ~3544000
################################################################################

for (i in seq(from = 3.2, to = 3.9, by = 0.1))
{
  n <- i * 10^6
  print(n)
  data <- build_data(n, nb_seg)
  g <- graph(type = "std", penalty = 2*log(n))
  print(system.time(gfpop(data, mygraph = g, type = "poisson")))
}

### response: ___


################################################################################
##### DUST -- extrapolated guess ~30890000
################################################################################

for (i in seq(from = 2.9, to = 3.3, by = 0.05))
{
  n <- i * 10^7
  print(n)
  data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")
  print(system.time(dust.1D(data, 2*log(n), model = "poisson", method = "DUST", backend = "highway")))
}

### response: ___
