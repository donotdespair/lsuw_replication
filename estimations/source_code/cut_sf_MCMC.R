
cut_sf_MCMC    = function(sf, credMass = .99) {
  # sf - SF posterior obtained by applying function TranstoMR
  S         = dim(sf)[1]
  N         = dim(sf)[2]
  within_s  = array(NA, c(N, S))
  range_    = apply(sf, 2, HDInterval::hdi, credMass = credMass)
  
  for (j in 1:N) {
    within_s[j,]  = sf[,j] > range_[1,j] & sf[,j] < range_[2,j]
  }
  
  within    = apply(within_s, 2, all)
  
  sf        = sf[within,]
  
  cat("remained: ", 100 * sum(within) / S)
  
  return(sf)
}
