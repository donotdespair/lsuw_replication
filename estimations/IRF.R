
# IRF: of gdp to unanticipated tax shock over samples
############################################

library(bsvars)

# load("spartan/results/tax06n.rda")
# irfs06 = compute_impulse_responses(post, horizon = 20)

# load("spartan/results/taxMRn.rda")
# irfsMR = compute_impulse_responses(post, horizon = 20)

# load("spartan/results/tax23.rda")
# irfs23 = compute_impulse_responses(post, horizon = 20)

data("us_fiscal_lsuw")

# rm("post","sddr", spec)

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
ylim = c(-4, 4)


# # responses of gdp to be discussed in the paper
# #######################################################
# apply(multiplier * (1/median(irfs23[1,1,1,])) * irfs23[3,1,,], 1, median)
# apply(multiplier * (1/median(irfs23[1,1,1,])) * irfs23[3,1,,], 1, HDInterval::hdi, credMass = 0.68)
# 
# # responses of gdp across three samples
# #######################################################
# pdf(
#   file = "spartan/IRF.pdf", 
#   height = 3, 
#   width = 12
# )
# par(
#   mfrow = c(1,3),
#   mar = c(4,5,1.5,0)
# )
# bsvarTVPs::ribbon_plot(
#   multiplier * (1/median(irfs23[1,1,1,])) * irfs23[3,1,,], 
#   probability = .68, 
#   col = cols[1],
#   ylim = ylim,
#   main = "2023 sample",
#   ylab = "gdp response [%]",
#   xlab = "propagation horizon [quarters]",
#   cex.lab = 1.6,
#   cex.axis = 1.6,
#   cex.main = 1.6,
#   bty = "n",
#   lwd = 2
# )
# abline(h = 0)
# 
# bsvarTVPs::ribbon_plot(
#   multiplier * (1/median(irfs06[1,1,1,])) * irfs06[3,1,,], 
#   probability = .68, 
#   col = cols[3],
#   ylim = ylim,
#   main = "2006 sample",
#   xlab = "propagation horizon [quarters]",
#   cex.lab = 1.6,
#   cex.axis = 1.6,
#   cex.main = 1.6,
#   bty = "n",
#   lwd = 2
# )
# abline(h = 0)
# 
# # MRdata
# bsvarTVPs::ribbon_plot(
#   multiplier * (1/median(irfsMR[1,1,1,])) * irfsMR[3,1,,], 
#   probability = .68, 
#   col = cols[5],
#   ylim = ylim,
#   main = "MR sample",
#   xlab = "propagation horizon [quarters]",
#   cex.lab = 1.6,
#   cex.axis = 1.6,
#   cex.main = 1.6,
#   bty = "n",
#   lwd = 2
# )
# abline(h = 0)
# dev.off()
# 
# system("cp spartan/IRF.pdf /Users/twozniak/Research/Carnaval-gulasz-Fei/fei-tomasz-sv-prior/cg_paper/graph_IRF.pdf")





# responses of gdp across three paper
#######################################################

bp_irf  = read.csv("spartan/BP_reproductions/IRF_BP1A.csv")
mr_irf  = read.csv("spartan/BP_reproductions/IRF_MR4A.csv")
le_irf  = read.csv("spartan/BP_reproductions/IRF_LE.csv")
hh      = nrow(mr_irf)

ylim    = c(-3, 5)
cols    = c("darkorchid1",
            adjustcolor("darkorchid1", alpha.f = 0.3), 
            "orchid1",
            adjustcolor("orchid1", alpha.f = 0.3), 
            "orchid4",
            adjustcolor("orchid4", alpha.f = 0.3), 
            "darkorchid4",
            adjustcolor("darkorchid4", alpha.f = 0.3)
)

load("spartan/results/tax06nPM.rda")
irfs06 = compute_impulse_responses(post, horizon = 20)


pdf(
  file = "spartan/IRF_comparizon.pdf",
  height = 9,
  width = 11
)
par(
  mfrow = c(2,2),
  mar = c(4,5,1.5,0)
)
bsvarTVPs::ribbon_plot(
  multiplier * (1/median(irfs06[1,1,1,])) * irfs06[3,1,,], 
  probability = .9, 
  col = cols[1],
  ylim = ylim,
  main = "2006 sample",
  ylab = "gdp response [%]",
  xlab = "",
  cex.lab = 1.6,
  cex.axis = 1.6,
  cex.main = 1.6,
  bty = "n",
  lwd = 2
)
abline(h = 0)

# BP results
plot(
  x = 0:20,
  y = rep(0, 21),
  main = "BP results",
  ylim = ylim,
  ylab = "",
  xlab = "",
  cex.lab = 1.6,
  cex.axis = 1.6,
  cex.main = 1.6,
  type = "n",
  bty = "n",
  axes = FALSE
)
polygon(
  x = c(1:hh, rev(1:hh)),
  y = c(bp_irf$lower, rev(bp_irf$upper)),
  col = cols[4],
  border = NA
)
lines(
  x = 1:hh,
  y = bp_irf$mean,
  col = cols[3],
  lwd = 2
)
abline(h = 0)
axis(1, c(0, 5, 10, 15, 20), c(0, 5, 10, 15, 20), cex.axis = 1.6)
axis(2, c(-2, 0, 2, 4), c(-2, 0, 2, 4), cex.axis = 1.6)

# MR results
plot(
  x = 0:20,
  y = rep(0, 21),
  main = "MR results",
  ylim = ylim,
  ylab = "gdp response [%]",
  xlab = "propagation horizon [quarters]",
  cex.lab = 1.6,
  cex.axis = 1.6,
  cex.main = 1.6,
  type = "n",
  bty = "n",
  axes = FALSE
)
polygon(
  x = c(1:hh, rev(1:hh)),
  y = c(mr_irf[,2], rev(mr_irf[,3])),
  col = cols[6],
  border = NA
)
lines(
  x = 1:hh,
  y = mr_irf[,1],
  col = cols[5],
  lwd = 2
)
abline(h = 0)
axis(1, c(0, 5, 10, 15, 20), c(0, 5, 10, 15, 20), cex.axis = 1.6)
axis(2, c(-2, 0, 2, 4), c(-2, 0, 2, 4), cex.axis = 1.6)

# LE results
plot(
  x = 0:20,
  y = rep(0, 21),
  main = "LE results",
  ylim = ylim,
  ylab = "",
  xlab = "propagation horizon [quarters]",
  cex.lab = 1.6,
  cex.axis = 1.6,
  cex.main = 1.6,
  type = "n",
  bty = "n",
  axes = FALSE
)
polygon(
  x = c(0:18, rev(0:18)),
  y = c(le_irf[,2], rev(le_irf[,3])),
  col = cols[8],
  border = NA
)
lines(
  x = 0:18,
  y = le_irf[,1],
  col = cols[7],
  lwd = 2
)
abline(h = 0)
axis(1, c(0, 5, 10, 15, 20), c(0, 5, 10, 15, 20), cex.axis = 1.6)
axis(2, c(-2, 0, 2, 4), c(-2, 0, 2, 4), cex.axis = 1.6)

dev.off()

# system("cp spartan/IRF_comparizon.pdf /Users/twozniak/Research/Carnaval-gulasz-Fei/fei-tomasz-sv-prior/cg_paper/graph_IRF_comparizon.pdf")




