
check_restriction = function(post) {
  S       = dim(post$posterior$sigma2_omega)[2]
  N       = dim(post$posterior$sigma2_omega)[1]
  check   = matrix(NA, S, N)
  for (s in 1:S) {
    check[s,] = post$posterior$sigma2_omega[,s] / (1 - post$posterior$rho[,s]^2) <= 1
  }
  out     = apply(
    cbind(apply(check,1,all), check), 
    2, mean
  )
  names(out) = c("joint", 1:N)
  return(out) 
}
