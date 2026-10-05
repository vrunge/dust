# Keep the former capitalized spelling as a compatibility alias.
.dust_backend <- function(backend) {
  if (identical(backend, "Highway")) backend <- "highway"
  match.arg(backend, c("highway", "scalar"))
}
