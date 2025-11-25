# Extracted from test-dust.1D.object.R:66

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "dust", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)
library(dust)

# test -------------------------------------------------------------------------
obj_dust <- dust.object.1D("poisson", "randIndex_Eval4", 9)
penalty <- 2
data_all <- NULL
for(i in 1:500)
            {
              data <- rpois(1, 100)
              data_all <- c(data_all, data)
              obj_dust$append(data, penalty)
              obj_dust$update_partition()
            }
