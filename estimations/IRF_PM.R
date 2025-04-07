
# IRF: of gdp to unanticipated tax shock over samples
############################################

library(bsvars)

# load("spartan/results/tax06n.rda")
# irfs06 = compute_impulse_responses(post, horizon = 20)
# 
# load("spartan/results/taxMRn.rda")
# irfsMR = compute_impulse_responses(post, horizon = 20)

load("spartan/results/tax23nPM.rda")
irfs23 = compute_impulse_responses(post, horizon = 20)

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
# multiplier = -(data2006q4[3]/data2006q4[1])
ylim = c(-4, 4)


# impulse responses
# 23 sample with PM normalisation
#######################################################
pdf(
  file = "spartan/IRF_PM.pdf", 
  height = 7, 
  width = 12
)
par(
  mfrow = c(3,3)
)

par(
  mar = c(4,4.5,1.5,0)
)
for (i in 1:3) {
  for (j in 1:3) {
    if (j == 1) {
      ylab  = paste(colnames(us_fiscal_lsuw)[i], "response")
    } else {
      ylab  = ""
    }
    if (i == 3 & j == 2) {
      xlab = "forecast horizon [quarters]"
    } else {
      xlab = ""
    }
    if (i == 1) {
      main = substitute(w[i], list(i = j))
    } else {
      main = ""
    }
    bsvarTVPs::ribbon_plot(
      irfs23[i,j,,], 
      probability = .68, 
      col = cols[1],
      main = main,
      ylab = ylab,
      xlab = xlab,
      cex.lab = 1.6,
      cex.axis = 1.6,
      cex.main = 1.6,
      bty = "n",
      lwd = 2
    )
    abline(h = 0)
  }
}

dev.off()

system("cp spartan/IRF_PM.pdf /Users/twozniak/Research/Carnaval-gulasz-Fei/fei-tomasz-sv-prior/cg_paper/graph_IRF_PM.pdf")
