# Extracted from test-dust.1D.object.R:42

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "dust", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(testthat)
library(dust)

# test -------------------------------------------------------------------------
obj_dust <- dust.object.1D("variance", "randIndex_Eval4", 5)
penalty <- 2*log(5)
data_all <- NULL
for(i in 1:5)
            {
              data <- dataGenerator_1D(chpts = c(500,1000), parameters = c(1,1.5), type = "variance")
              data_all <- c(data_all, data)
              obj_dust$append(data, penalty)
              obj_dust$update_partition()
            }
