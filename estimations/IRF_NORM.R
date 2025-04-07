
# IRF: of gdp to unanticipated tax shock over samples
############################################

library(bsvars)

# load("spartan/results/tax06n.rda")
# irfs06 = compute_impulse_responses(post, horizon = 20)
# 
# load("spartan/results/taxMRn.rda")
# irfsMR = compute_impulse_responses(post, horizon = 20)
# 
# load("spartan/results/tax23.rda")
# irfs23 = compute_impulse_responses(post, horizon = 20)
# 
# 
# 
# load("spartan/results/tax06nMR.rda")
# irfs06mr = compute_impulse_responses(post, horizon = 20)
# 
# load("spartan/results/taxMRnMR.rda")
# irfsMRmr = compute_impulse_responses(post, horizon = 20)
# 
# load("spartan/results/tax23nMR.rda")
# irfs23mr = compute_impulse_responses(post, horizon = 20)



load("spartan/results/tax06nPM.rda")
irfs06pm = compute_impulse_responses(post, horizon = 20)

# load("spartan/results/taxMRnPM.rda")
# irfsMRpm = compute_impulse_responses(post, horizon = 20)

load("spartan/results/tax23nPM.rda")
irfs23pm = compute_impulse_responses(post, horizon = 20)




data("us_fiscal_lsuw")

rm("post","sddr", spec)

# colors
cols    = c("darkorchid1",
            adjustcolor("darkorchid1", alpha.f = 0.3), 
            "orchid4",
            adjustcolor("orchid4", alpha.f = 0.3), 
            "darkorchid4",
            adjustcolor("darkorchid4", alpha.f = 0.3)
)


data2006q4 = us_fiscal_lsuw[which(zoo::index(us_fiscal_lsuw) == 2006.75), ]
multiplier = -5 * (data2006q4[3]/data2006q4[1])
ylim = c(-2, 2)

# computations for interpretations
#######################################################
apply(multiplier * (1/median(irfs23pm[1,1,1,])) * irfs23pm[3,1,,], 1, median)
apply(multiplier * (1/median(irfs06pm[1,1,1,])) * irfs06pm[3,1,,], 1, median)

# responses of gdp
#######################################################
pdf(
  file = "spartan/IRF_NORM.pdf",
  height = 5,
  width = 11
)
par(
  mfrow = c(1,2)
)

# 23 sample
par(
  mar = c(4,5,1.5,0)
)

# PM-ordering

bsvarTVPs::ribbon_plot(
  multiplier * (1/median(irfs23pm[1,1,1,])) * irfs23pm[3,1,,], 
  probability = .68, 
  col = cols[1],
  ylim = ylim,
  main = "2023 sample",
  ylab = "gdp response [%]",
  xlab = "propagation horizon [quarters]",
  cex.lab = 1.6,
  cex.axis = 1.6,
  cex.main = 1.6,
  bty = "n",
  lwd = 2
)
abline(h = 0)

bsvarTVPs::ribbon_plot(
  multiplier * (1/median(irfs06pm[1,1,1,])) * irfs06pm[3,1,,],
  probability = .68,
  col = cols[3],
  ylim = ylim,
  main = "2006 sample",
  xlab = "propagation horizon [quarters]",
  cex.lab = 1.6,
  cex.axis = 1.6,
  cex.main = 1.6,
  bty = "n",
  lwd = 2
)
abline(h = 0)

dev.off()

# system("cp spartan/IRF_NORM.pdf /Users/twozniak/Research/Carnaval-gulasz-Fei/fei-tomasz-sv-prior/cg_paper/graph_IRF_NORM.pdf")
