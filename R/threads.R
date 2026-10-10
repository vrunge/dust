## all the cores for the methods keeping many indices (2 under R CMD check --as-cran)
.default_threads <- function(method)
{
  if (!method[1] %in% c("OP", "PELT", "PELTpar")) return(1L)
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) return(2L)
  max(1L, parallel::detectCores(), na.rm = TRUE)
}

## NA: no size (binom and negbin data already divided by it)
.size_arg <- function(size)
{
  if (is.null(size)) return(NA_real_)
  if (!is.numeric(size) || length(size) != 1L || !is.finite(size) || size <= 0)
    stop("size must be a positive number", call. = FALSE)
  as.double(size)
}
