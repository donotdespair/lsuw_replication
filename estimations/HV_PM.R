
# HV: heteroskedasticity verification plots
############################################

library(bsvars)

load("spartan/results/tax23nPM.rda")
post23 = post

load("spartan/results/tax06nPM.rda")
post06 = post

load("spartan/results/taxMRnPM.rda")
postMR = post

data("us_fiscal_lsuw")

rm("post","sddr","spec")

# colors
cols    = c("darkorchid1",
            adjustcolor("darkorchid1", alpha.f = 0.3), 
            "orchid4",
            adjustcolor("orchid4", alpha.f = 0.3), 
            "darkorchid4",
            adjustcolor("darkorchid4", alpha.f = 0.3)
          )

# # PRELIMINARY: conditional sds
# csds23    = compute_conditional_sd(post23)
# csds06    = compute_conditional_sd(post06)
# 
# par(mfrow = c(3,1))
# for (j in 1:3) {
#   bsvarTVPs::ribbon_plot(csds06[j,,]); abline(h = 1, col = 1)
# }

# conditional vars
dates23       = zoo::index(us_fiscal_lsuw)[-(1:4)]
dates06       = dates23[-((length(dates23) - 66):length(dates23))]
datesMR       = dates06[-(1:8)]

sigma223      = post23$posterior$sigma^2
sigma223_med  = apply(sigma223, 1:2, mean)
sigma223_hdi  = apply(sigma223, 1:2, HDInterval::hdi, credMass = .90)

sigma206      = post06$posterior$sigma^2
sigma206_med  = apply(sigma206, 1:2, mean)
sigma206_hdi  = apply(sigma206, 1:2, HDInterval::hdi, credMass = .90)

sigma2MR      = postMR$posterior$sigma^2
sigma2MR_med  = apply(sigma2MR, 1:2, mean)
sigma2MR_hdi  = apply(sigma2MR, 1:2, HDInterval::hdi, credMass = .90)

# investigating the dates when variances are significantly greater than 1
############################################

# 2023 sample
dates23[which(sigma223_hdi[1,1,]>1)] # ttr
dates23[which(sigma223_hdi[1,2,]>1)] # gs
dates23[which(sigma223_hdi[1,3,]>1)] # gdp

# 2006 sample
dates06[which(sigma206_hdi[1,1,]>1)] # ttr
dates06[which(sigma206_hdi[1,2,]>1)] # gs
dates06[which(sigma206_hdi[1,3,]>1)] # gdp

# MR data
datesMR[which(sigma2MR_hdi[1,1,]>1)] # ttr
datesMR[which(sigma2MR_hdi[1,2,]>1)] # gs
datesMR[which(sigma2MR_hdi[1,3,]>1)] # gdp


# the plot
############################################
plot_range    = matrix(NA, 3, 2)
for (i in 1:3) {
  plot_range[i,]  = range(c(sigma223_hdi[,i,], sigma206_hdi[,i,], sigma2MR_hdi[,i,]))
}

pdf(
  file = "spartan/HV_PM.pdf",
  height  = 2.5,
  width   = 12
)
par(
  mfcol = c(1,3)
)
for (n in 1:3) {
  
  # 23 plots
  par(mar = c(4,5,1,0))
  ylab = ""
  if (n == 1) {
    ylab = paste("2023 sample")
    par(mar = c(4,5,1,0))
  }
  plot(
    x = dates23,
    y = sigma223_med[n,],
    type = "l",
    ylim = range(sigma223_hdi[,n,]),
    xlab = "",
    ylab = ylab,
    main = paste("conditional variance of ", colnames(us_fiscal_lsuw)[n]),
    col = "white",
    bty = "n",
    cex.lab = 1.6,
    cex.axis = 1.6,
    cex.main = 1.6,
  )
  polygon(
    x = c(dates23, rev(dates23)),
    y = c(sigma223_hdi[1,n,], rev(sigma223_hdi[2,n,])),
    col = cols[2],
    border = NA
  )
  lines(
    x = dates23,
    y = sigma223_med[n,],
    lwd = 2,
    col = cols[1]
  )
  abline(h = 1, col = 1)
  # 
  # # 06 plots
  # par(mar = c(4,5,1,0))
  # ylab = ""
  # if (n == 1) {
  #   ylab = paste("2006 sample")
  #   par(mar = c(4,5,1,1))
  # }
  # plot(
  #   x = dates23,
  #   y = sigma223_med[n,],
  #   type = "l",
  #   ylim = range(sigma206_hdi[,n,]),
  #   xlab = "",
  #   ylab = ylab,
  #   col = "white",
  #   bty = "n",
  #   cex.lab = 1.6,
  #   cex.axis = 1.6,
  #   cex.main = 1.6,
  # )
  # polygon(
  #   x = c(dates06, rev(dates06)),
  #   y = c(sigma206_hdi[1,n,], rev(sigma206_hdi[2,n,])),
  #   col = cols[4],
  #   border = NA
  # )
  # lines(
  #   x = dates06,
  #   y = sigma206_med[n,],
  #   lwd = 2,
  #   col = cols[3]
  # )
  # abline(h = 1, col = 1)
  # 
  # # MR plots
  # par(mar = c(4,5,1,0))
  # ylab = ""
  # if (n == 1) {
  #   ylab = paste("MR sample")
  #   par(mar = c(4,5,1,0))
  # }
  # plot(
  #   x = dates23,
  #   y = sigma223_med[n,],
  #   type = "l",
  #   ylim = range(sigma2MR_hdi[,n,]),
  #   xlab = "time",
  #   ylab = ylab,
  #   col = "white",
  #   bty = "n",
  #   cex.lab = 1.6,
  #   cex.axis = 1.6,
  #   cex.main = 1.6,
  # )
  # polygon(
  #   x = c(datesMR, rev(datesMR)),
  #   y = c(sigma2MR_hdi[1,n,], rev(sigma2MR_hdi[2,n,])),
  #   col = cols[4],
  #   border = NA
  # )
  # lines(
  #   x = datesMR,
  #   y = sigma2MR_med[n,],
  #   lwd = 2,
  #   col = cols[3]
  # )
  # abline(h = 1, col = 1)
}
dev.off()

system("cp spartan/HV_PM.pdf /Users/twozniak/Research/Carnaval-gulasz-Fei/fei-tomasz-sv-prior/cg_paper/graph_HV_PM.pdf")
