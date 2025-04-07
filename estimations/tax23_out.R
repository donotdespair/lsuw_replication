
library(bsvars)
source("spartan/source_code/TransformtoMR.R")
source("spartan/source_code/check_restriction.R")
Rcpp::sourceCpp("spartan/source_code/normalise_jaro.cpp")

# Post processing
# This model works without normalisation
########################################################
data("us_fiscal_lsuw")
plot(us_fiscal_lsuw)

load("spartan/results/tax23.rda")

sddr$logSDDR
sddr$log_SDDR_se
BF = 1/exp(sddr$logSDDR); BF
100 * BF / (1 + BF) # posterior probability of heteroskedasticity

check_restriction(post)

plot.ts(t(post$posterior$B[1,,]))
plot.ts(t(post$posterior$B[2,,]))
plot.ts(t(post$posterior$B[3,,]))

apply(post$posterior$B, 1:2, mean)
apply(post$posterior$B, 1:2, sd)

apply(post$posterior$A, 1:2, mean)
apply(post$posterior$A, 1:2, sd)

ss <- post |> compute_structural_shocks()
par(mfrow = c(3,1))
bsvarTVPs::ribbon_plot(ss[1,,])
bsvarTVPs::ribbon_plot(ss[2,,])
bsvarTVPs::ribbon_plot(ss[3,,])

irfs <- post |> compute_impulse_responses(horizon = 20)
par(mfrow = c(3,3))
for (i in 1:3) {
  for (j in 1:3) {
    bsvarTVPs::ribbon_plot(irfs[i,j,,]); abline(h = 0, col = 1)
  }
}

csds <- post |> compute_conditional_sd()
par(mfrow = c(3,1))
for (j in 1:3) {
  bsvarTVPs::ribbon_plot(csds[j,,]); abline(h = 1, col = 1)
}

post_sf   = TranstoMR(post$posterior)
apply(post_sf, 2, mean)
apply(post_sf, 2, sd)
apply(post_sf, 2, HDInterval::hdi, credMass = .68)



# normalization a la Lewis 2021
########################################################
Rcpp::sourceCpp("spartan/source_code/normalise_jaro.cpp")
sample = "23"

# MR results for Proxy VAR
ThetaG = -0.2
ThetaY = 3.13
SigmaT = 2.54/100
GammaT = 0.06
GammaY = 0
SigmaG = 2.35/100
ZetaT  = -0.36
ZetaG  = 0.1
SigmaY = 1.54/100

b11 = 1/(SigmaT*(1 - ThetaG*GammaT))
b12 = -ThetaG/(SigmaT*(1 - ThetaG*GammaT))
b13 = (GammaY*ThetaG - ThetaY)/(SigmaT*(1 - ThetaG*GammaT))
b21 = -GammaT/(SigmaG*(1 - ThetaG*GammaT))
b22 = 1/(SigmaG*(1 - ThetaG*GammaT))
b23 = (GammaT*ThetaY - GammaY)/(SigmaG*(1 - ThetaG*GammaT))
b31 = -ZetaT/SigmaY
b32 = -ZetaG/SigmaY
b33 = 1/SigmaY
B0hat_MR    = matrix(c(b11,b21,b31,b12,b22,b32,b13,b23,b33),3,3)

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
B0hat_BP    = matrix(c(b11,b21,b31,b12,b22,b32,b13,b23,b33),3,3)


load(paste0("spartan/results/tax", sample,".rda"))
post_normJ_BP     = normalise_jaro_bsvar_sv(post$posterior, B0hat_BP)
post              = bsvars::specify_posterior_bsvar_sv$new(post$last_draw, post_normJ_BP)
post$set_normalised()
sddr              = bsvars::verify_volatility(post)
save(post, spec, sddr, B0hat_BP, file = paste0("spartan/results/tax", sample,"nBP.rda"))

load(paste0("spartan/results/tax", sample,".rda"))
qq_normL_MR       = normalise_jaro_bsvar_sv(post$posterior, B0hat_MR)
post              = bsvars::specify_posterior_bsvar_sv$new(post$last_draw, qq_normL_MR)
post$set_normalised()
sddr              = bsvars::verify_volatility(post)
save(post, spec, sddr, B0hat_MR, file = paste0("spartan/results/tax", sample,"nMR.rda"))

load(paste0("spartan/results/tax", sample,".rda"))
B0hat_PM    = post$last_draw$starting_values$B
qq_normL_PM       = normalise_jaro_bsvar_sv(post$posterior, B0hat_PM)
post              = bsvars::specify_posterior_bsvar_sv$new(post$last_draw, qq_normL_PM)
post$set_normalised()
sddr              = bsvars::verify_volatility(post)
save(post, spec, sddr, B0hat_PM, file = paste0("spartan/results/tax", sample,"nPM.rda"))
