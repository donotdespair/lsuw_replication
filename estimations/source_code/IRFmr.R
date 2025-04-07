IRF <-function(qq,p=4,H=20,VAR.no=3,shock.no=1){
n  = nrow(qq$B)
S = dim(qq$B)[3]
IM        <- data.frame(cbind(seq(1,H)))
# IM        <- data.frame(cbind(seq(0,H)))
names(IM) <- "quarters"

qqpermute = qq


qqpermutesf <- TranstoMR(qqpermute)

# The issue is from B0[3,1]，because B0[3,1]=-ZetaT/SigmaY，ZT is positive and significant here,but negative in MR BP

  # IRF = array(NA,dim = c(3,H+1,S))
IRF = array(NA,dim = c(3,H,S))
  DAT.TRY = 0.1746
  VAR.tshocksize =  0.01/DAT.TRY
  for (x in 1:S) {
    ZetaT   = array(NA, c(1,S))
    ZetaG   = array(NA, c(1,S))
    ThetaY = array(NA, c(1,S))
    ThetaG = array(NA, c(1,S))
    GammaY = array(NA, c(1,S))
    GammaT = array(NA, c(1,S))
    SigmaT = array(NA, c(1,S))
    SigmaG = array(NA, c(1,S))
    SigmaY = array(NA, c(1,S))
    Amr = cbind(c(1,0,-qqpermutesf$ZetaT[x]),c(0,1,-qqpermutesf$ZetaG[x]),c(-qqpermutesf$ThetaY[x],-qqpermutesf$GammaY[x],1))
    Bmr = diag(3)
    Bmr[2,1] = qqpermutesf$GammaT[x]
    Bmr[1,2] = qqpermutesf$ThetaG[x]
    Dmr      = diag(c(qqpermutesf$SigmaT[x],qqpermutesf$SigmaG[x],qqpermutesf$SigmaY[x]))
    #A0 computed this way is the same as normal way
    A0       = solve(Amr)%*%Bmr%*%Dmr
    A   = qqpermute$A[,,x]
    # A0  = solve(qqpermute$B[,,x])
    F                        = matrix(0,n*p,n*p)
    F[(n+1):(n*p),1:(n*p-n)] = diag(n*p-n)
    F[1:n,1:(n*p)]           = A[,1:(n*p)]
    #IRF
    # yy             = array(NA,dim = c(3,H+1))
    yy             = array(NA,dim = c(3,H))
    J     = cbind(diag(n), matrix(0,n,n*(p-1)))
    Ftemp = F%^%0
    # Ftemp = F%^%1
    # for (ii in 1:(H+1)) {
    for (ii in 1:(H)) {
      ytemp = J%*%(Ftemp)%*%t(J)%*%A0[,shock.no]
      yy[,ii] = -100*ytemp/A0[1,1]*VAR.tshocksize
      Ftemp = Ftemp%*%F
    }
    #standardize the shock 
    IRF[,,x] = yy
  }
  Irfall     <- apply(IRF,1:2,median)
  Irfbandall <- apply(IRF,1:2,HDInterval::hdi,credMass = 0.95)
  
  Irf.tax1     <- Irfall[VAR.no,]
  Irfband.tax1 <- Irfbandall[,VAR.no,] 
  
  IRF          <- data.frame(Irf.tax1,t(Irfband.tax1))
  IM           <- cbind(IM,IRF)
  return(IM)
}