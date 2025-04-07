
normalizationL2021 <-function(posterior,B0hatBP){
  N = nrow(posterior$h)
  S = ncol(posterior$rho)
  B0hat <- posterior$B[,,S]
  # B0hat <- diag(sign(diag(B0hat))) %*% B0hat
  IRFhat <- solve(B0hat)%*%diag(1/diag(solve(B0hat)))
  A0hat  <- solve(B0hat)
  posterior$B0 = array(NA,dim = c(N,N,S))
  posterior$A0 = array(NA,dim = c(N,N,S))
  posterior$restriction24 = array(NA,dim = c(N,S))
###############################################################################  
  pcM        <- combinat::permn(N)
  P          <-lapply(1:factorial(N), function(n){diag(N)[pcM[[n]],]})
  Sdiag      <-expand.grid(
    c(1,-1),
    c(1,-1),
    c(1,-1)
  )
    for (s in 1:S) {
      B0s      <- posterior$B[,,s]
    allB0      <- list()
    allA0      <- list()
    #PS is the permutation of rows and signs
    allPS      <- list() 
    metrics    <-vector()
      for (n in 1:2^N) {
        Smatrix    <-diag(Sdiag[n,])
        B0n        <-lapply(1:factorial(N), function(x){diag(N)[pcM[[x]],]%*%Smatrix%*%B0s})
        scores     <-unlist(lapply(1:(factorial(N)), function(x){t(as.vector(B0n[[x]]-B0hat))%*%as.vector(B0n[[x]]-B0hat)}))
        
        #############################################################################################################
        A0n         <-lapply(1:factorial(N), function(x){solve(B0n[[x]])})
        # scores     <-unlist(lapply(1:(factorial(N)), function(x){t(as.vector(A0n[[x]]-A0hat))%*%as.vector(A0n[[x]]-A0hat)}))
        ################################################################################################################
        # IRFn.at0    <-lapply(1:factorial(N), function(x){solve(B0n[[x]])%*%diag(1/diag(solve(B0n[[x]])))})
        # scores      <-unlist(lapply(1:(factorial(N)), function(x){t(as.vector(IRFn.at0[[x]]-IRFhat))%*%as.vector(IRFn.at0[[x]]-IRFhat)}))
        # ###############################################################################################
        B0n.lowest <-B0n[which.min(scores)]
        A0n.lowest <-A0n[which.min(scores)]
        PS.lowest  <-list((diag(N)[pcM[[which.min(scores)]],]%*%Smatrix))
        allB0      <-c(allB0,B0n.lowest)
        allA0      <-c(allA0,A0n.lowest)
        allPS      <-c(allPS,PS.lowest)
        metrics    <-c(metrics,scores[which.min(scores)])
      }

    #update MCMC
    PS         <- allPS[[which.min(metrics)]] 
    B0s        <- allB0[[which.min(metrics)]]
    A0s        <- allA0[[which.min(metrics)]]
    # posterior$Bplus[,,s]       <- PS%*%posterior$Bplus[,,s]
    # posterior$A[,,s]       <- PS%*%posterior$A[,,s]
    posterior$B0[,,s]          <- B0s
    posterior$A0[,,s]          <- A0s
    posterior$h[,,s]           <- abs(PS)%*%posterior$h[,,s]
    posterior$sigma[,,s]       <- abs(PS)%*%posterior$sigma[,,s]        
    posterior$rho[,s]          <- abs(PS)%*%posterior$rho[,s]
    posterior$omega[,s]        <- abs(PS)%*%posterior$omega[,s]
    posterior$sigma2_omega[,s] <- abs(PS)%*%posterior$sigma2_omega[,s]
    posterior$S[,,s]           <- abs(PS)%*%posterior$S[,,s]
    posterior$restriction24[,s]      <- posterior$sigma2_omega[,s]/(1-posterior$rho[,s]^2)<=1
    }
  
  ##########################################################################################
  B0mean <- apply(posterior$B0,1:2,mean)
  B0hat   <- B0hatBP
  ##########################################################################################
  
  pcM        <- combinat::permn(N)
  P          <-lapply(1:factorial(N), function(n){diag(N)[pcM[[n]],]})
  Sdiag      <-expand.grid(
    c(1,-1),
    c(1,-1),
    c(1,-1)
  )
  
    # B0s <- posterior$B0[,,s]
    B0s        <- B0mean
    allB0      <- list()
    allA0      <- list()
    #PS is the permutation of rows and signs
    allPS      <- list() 
    metrics    <-vector()
    for (n in 1:2^N) {
      Smatrix    <-diag(Sdiag[n,])
      B0n        <-lapply(1:factorial(N), function(x){diag(N)[pcM[[x]],]%*%Smatrix%*%B0s})
      scores     <-unlist(lapply(1:(factorial(N)), function(x){t(as.vector(B0n[[x]]-B0hat))%*%as.vector(B0n[[x]]-B0hat)}))
      
      #############################################################################################################
      # A0n         <-lapply(1:factorial(N), function(x){solve(B0n[[x]])})
      # scores      <-unlist(lapply(1:(factorial(N)), function(x){t(as.vector(A0n[[x]]-A0hat))%*%as.vector(A0n[[x]]-A0hat)}))
      ################################################################################################################
      # IRFn.at0    <-lapply(1:factorial(N), function(x){solve(B0n[[x]])%*%diag(1/diag(solve(B0n[[x]])))})
      # scores      <-unlist(lapply(1:(factorial(N)), function(x){t(as.vector(IRFn.at0[[x]]-IRFhat))%*%as.vector(IRFn.at0[[x]]-IRFhat)}))
      # #############################################################################################################
      B0n.lowest <-B0n[which.min(scores)]
      A0n.lowest <-A0n[which.min(scores)]
      PS.lowest  <-list((diag(N)[pcM[[which.min(scores)]],]%*%Smatrix))
      allB0      <-c(allB0,B0n.lowest)
      allA0      <-c(allA0,A0n.lowest)
      allPS      <-c(allPS,PS.lowest)
      metrics    <-c(metrics,scores[which.min(scores)])
    }
    
    #update MCMC
    PS         <- allPS[[which.min(metrics)]] 
    B0s        <- allB0[[which.min(metrics)]]
    A0s        <- allA0[[which.min(metrics)]]
    for (s in 1:S) {
      posterior$B[,,s]           <- PS%*%posterior$B0[,,s]
      posterior$B0[,,s]          <- PS%*%posterior$B0[,,s]
      posterior$A0[,,s]          <- solve(posterior$B0[,,s])
      posterior$h[,,s]           <- abs(PS)%*%posterior$h[,,s]
      posterior$sigma[,,s]       <- abs(PS)%*%posterior$sigma[,,s]        
      posterior$rho[,s]          <- abs(PS)%*%posterior$rho[,s]
      posterior$omega[,s]        <- abs(PS)%*%posterior$omega[,s]
      posterior$sigma2_omega[,s] <- abs(PS)%*%posterior$sigma2_omega[,s]
      posterior$S[,,s]           <- abs(PS)%*%posterior$S[,,s]
      posterior$restriction24[,s]<- posterior$sigma2_omega[,s]/(1-posterior$rho[,s]^2)<=1
    }

  
  return(posterior)
}

