
thin_bsvar_sv <- function(qq, by = 100){
  T = dim(qq$posterior$S)[3]
  post = seq(1, T, by = by) 
  qq$posterior$B            = qq$posterior$B[,,post] 
  qq$posterior$A            = qq$posterior$A[,,post] 
  qq$posterior$hyper        = qq$posterior$hyper[,post]
  qq$posterior$h            = qq$posterior$h[,,post]
  qq$posterior$rho          = qq$posterior$rho[,post]
  qq$posterior$omega        = qq$posterior$omega[,post]
  qq$posterior$S            = qq$posterior$S[,,post]
  qq$posterior$sigma2_omega = qq$posterior$sigma2_omega[,post]
  qq$posterior$sigma        = qq$posterior$sigma[,,post]
  qq$posterior$s_           = qq$posterior$s_[,post]
  return(qq)
}

thin_sf <- function(qq, by = 100){
  T = dim(qq$S)[3]
  post = seq(1, T, by = by) 
  qq$B            = qq$B[,,post] 
  qq$A            = qq$A[,,post] 
  qq$hyper        = qq$hyper[,post]
  qq$h            = qq$h[,,post]
  qq$rho          = qq$rho[,post]
  qq$omega        = qq$omega[,post]
  qq$S            = qq$S[,,post]
  qq$sigma2_omega = qq$sigma2_omega[,post]
  qq$sigma        = qq$sigma[,,post]
  qq$s_           = qq$s_[,post]
  return(qq)
}

