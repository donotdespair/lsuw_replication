
cut_MCMC    = function(posterior, credMass = .99) {
  S         = dim(posterior$B)[3]
  N         = dim(posterior$B)[1]
  within_s  = array(NA, c(N,N, S))
  range_    = apply(posterior$B, 1:2, HDInterval::hdi, credMass = credMass)
  
  for (i in 1:N) {
    for (j in 1:N) {
    within_s[i,j,]  = posterior$B[i,j,] > range_[1,i,j] & posterior$B[i,j,] < range_[2,i,j]
    }
  }
  within    = apply(within_s, 3, all)
  
  posterior$B       = posterior$B[,,within]
  posterior$A       = posterior$A[,,within]
  posterior$hyper   = posterior$hyper[,within]
  posterior$h       = posterior$h[,,within]
  posterior$rho     = posterior$rho[,within]
  posterior$omega   = posterior$omega[,within]
  posterior$S       = posterior$S[,,within]
  posterior$sigma2_omega  = posterior$sigma2_omega[,within]
  posterior$s_      = posterior$s_[,within]
  posterior$sigma   = posterior$sigma[,,within]
  
  cat("remained: ", 100 * sum(within) / S)
  
  return(posterior)
}
