
library(bsvars)
source("spartan/source_code/TransformtoMR.R")
source("spartan/source_code/check_restriction.R")
Rcpp::sourceCpp("spartan/source_code/normalise_jaro.cpp")

# Post processing
# This model works without normalisation
########################################################
rm(list = ls())

# data("us_fiscal_lsuw")
# plot(us_fiscal_lsuw)

load("spartan/results/tax23ms06.rda")


plot.ts(t(apply(qqq$posterior$xi, 1:2, mean)))
plot.ts(t(qqq$posterior$xi[,,59506]))

apply(qqq$posterior$xi, 1, mean)

plot.ts(t(bbb[3,,1,]))
plot.ts(t(bbb[3,,2,]))

par(mfcol = c(3,2))
for ( n in 1:3) {
  for (m in 1:2) {
    hist(qqq$posterior$omega[3,2,], breaks = 200)
  }
}

bbb[,,,50000]










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
