# "Highway" is also accepted
.dust_backend <- function(backend) {
  if (identical(backend, "Highway")) backend <- "highway"
  match.arg(backend, c("highway", "scalar"))
}
