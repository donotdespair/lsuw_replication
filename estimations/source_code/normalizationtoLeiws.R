
normalizationLewis <-function(posterior,Hhat){
  N = nrow(posterior$h)
  S = ncol(posterior$rho)
  pcM        <- combinat::permn(N)
  P          <-lapply(1:factorial(N), function(n){diag(N)[pcM[[n]],]})
  posterior$SS <- array(NA,dim = c(1,s))
  Sdiag      <-expand.grid(
    c(1,-1),
    c(1,-1),
    c(1,-1)
  )
    for (s in 1:S) {
    B0s <- posterior$B0[,,s]
    allB0      <- list() 
    allPS      <- list() 
    metrics    <-vector()
    # permutes   <-lapply(1:factorial(N), function(n){diag(N)[pcM[[n]],]%*%B0s})
    # reverse    <-lapply(1:factorial(N), function(n){permutes[[n]]%*%(-diag(3))})
      for (n in 1:2^N) {
        Smatrix    <-diag(Sdiag[n,])
        B0n        <-lapply(1:factorial(N), function(x){diag(N)[pcM[[x]],]%*%Smatrix%*%B0s})
        H0         <-lapply(1:factorial(N), function(x){solve(B0n[[x]])%*%diag(1/diag(solve(B0n[[x]])))})
        scores     <-unlist(lapply(1:(factorial(N)), function(x){t(as.vector(H0[[x]]-Hhat))%*%as.vector(H0[[x]]-Hhat)}))
        B0n.lowest <-B0n[which.min(scores)]
        PS.lowest  <-list((diag(N)[pcM[[which.min(scores)]],]%*%Smatrix))
        allB0      <-c(allB0,B0n.lowest)
        allPS      <-c(allPS,PS.lowest)
        metrics    <-c(metrics,scores[which.min(scores)])
      }

    #update MCMC
    PS         <- allPS[[which.min(metrics)]] 
    B0s        <- allB0[[which.min(metrics)]]
    posterior$Bplus[,,s]       <- B0s%*%solve(posterior$B0[,,s])%*%posterior$Bplus[,,s]
    posterior$B0[,,s]          <- B0s
    posterior$h[,,s]           <- abs(PS)%*%posterior$h[,,s]
    posterior$sigma[,,s]       <- abs(PS)%*%posterior$sigma[,,s]        
    posterior$rho[,s]          <- abs(PS)%*%posterior$rho[,s]
    posterior$omega[,s]        <- abs(PS)%*%posterior$omega[,s]
    posterior$sigma2_omega[,s] <- abs(PS)%*%posterior$sigma2_omega[,s]
    posterior$S[,,s]           <- abs(PS)%*%posterior$S[,,s]
    posterior$SS[,s]           <- min(metrics)
    }
  
  return(posterior)
}

