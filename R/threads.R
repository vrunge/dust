## default number of threads: all the cores for the methods keeping many
## indices (OP, PELT, PELTpar), 2 under R CMD check --as-cran, 1 otherwise
.default_threads <- function(method)
{
  if (!method[1] %in% c("OP", "PELT", "PELTpar")) return(1L)
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) return(2L)
  max(1L, parallel::detectCores(), na.rm = TRUE)
}
