# as tax23.R but change prior for A from values c(10,10,10,10) to c(10,10,20,10)

model = "_10"

Rcpp::sourceCpp("source_code/normalise_jaro.cpp")

# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

library(bsvars)
load("MR2006.RData")

yMR   = ts(
          rbind(
            t(matrix(t(X)[1,1:12], 3,4)),
            t(Y)
          ), start = c(1950,1), frequency = 4
        )

tr = -(115.5:112.5)
exMR  = ts(
          rbind(
            cbind(tr, tr^2, rep(0,length(tr))),
            t(X)[,14:16]
          ), start = c(1950,1), frequency = 4
        )

set.seed(1234)
B         = matrix(TRUE, 3, 3)
spec      = specify_bsvar_sv$new(
              data = yMR,
              p    = 4,
              exogenous = exMR,
              B    = B
            )

spec$prior$hyper_s_AA = 20

burn      = estimate(spec, 3e5, thin = 1e4)
post      = estimate(burn, 6e5, thin = 10)

sddr      = verify_volatility(post)



# BP results for ThetaY = 2.08
ThetaG = -0.06
ThetaY = 2.08
SigmaT = 2.24/100
GammaT = 0
GammaY = 0
SigmaG = 2.36/100
ZetaT  = -0.08
ZetaG  = 0.07
SigmaY = 0.97/100

b11 = 1/(SigmaT*(1 - ThetaG*GammaT))
b12 = -ThetaG/(SigmaT*(1 - ThetaG*GammaT))
b13 = (GammaY*ThetaG - ThetaY)/(SigmaT*(1 - ThetaG*GammaT))
b21 = -GammaT/(SigmaG*(1 - ThetaG*GammaT))
b22 = 1/(SigmaG*(1 - ThetaG*GammaT))
b23 = (GammaT*ThetaY - GammaY)/(SigmaG*(1 - ThetaG*GammaT))
b31 = -ZetaT/SigmaY
b32 = -ZetaG/SigmaY
b33 = 1/SigmaY
B0hat_BP   = matrix(c(b11,b21,b31,b12,b22,b32,b13,b23,b33),3,3)

# compute the normalised posterior
########################################################
post_normJ_BP     = normalise_jaro_bsvar_sv(post$posterior, B0hat_BP)
post              = specify_posterior_bsvar_sv$new(post$last_draw, post_normJ_BP)
post$set_normalised()
sddr              = verify_volatility(post)
save(spec, post, sddr, file = paste0("results/taxMR", model,".rda"))
