TranstoMR <- function(qq) {
  S               = ncol(qq$rho)
  
  posteriors = data.frame(
    ZetaT    = rep(NA, S),
    ZetaG    = rep(NA, S),
    ThetaY   = rep(NA, S),
    ThetaG   = rep(NA, S),
    GammaY   = rep(NA, S),
    GammaT   = rep(NA, S),
    SigmaT   = rep(NA, S),
    SigmaG   = rep(NA, S),
    SigmaY   = rep(NA, S)
  )
  for (s in 1:S) {
    B0s     <- qq$B[,,s]
    ItaT    <- -B0s[3,1]/B0s[3,3]
    ItaG    <- -B0s[3,2]/B0s[3,3]
    ThetaY  <-(B0s[1,3]*B0s[2,2]-B0s[1,2]*B0s[2,3])/(B0s[1,2]*B0s[2,1]-B0s[1,1]*B0s[2,2])
    ThetaG  <- -B0s[1,2]/B0s[1,1]
    GammaY  <- (-B0s[1,3]*B0s[2,1]+B0s[1,1]*B0s[2,3])/(B0s[1,2]*B0s[2,1]-B0s[1,1]*B0s[2,2])
    GammaT  <- -B0s[2,1]/B0s[2,2]
    SigmaT  <- B0s[2,2]/(-B0s[1,2]*B0s[2,1]+B0s[1,1]*B0s[2,2])
    SigmaG  <- B0s[1,1]/(-B0s[1,2]*B0s[2,1]+B0s[1,1]*B0s[2,2])
    SigmaY  <- 1/B0s[3,3]
    
    posteriors$ZetaT[s]  <- ItaT
    posteriors$ZetaG[s]  <- ItaG
    posteriors$ThetaY[s] <- ThetaY
    posteriors$ThetaG[s] <- ThetaG
    posteriors$GammaY[s] <- GammaY
    posteriors$GammaT[s] <- GammaT
    posteriors$SigmaT[s] <- SigmaT
    posteriors$SigmaG[s] <- SigmaG
    posteriors$SigmaY[s] <- SigmaY
  }
  
  return(posteriors)
}

