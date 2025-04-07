
# OMEGA: omega posteriors for ttr with priors all samples
############################################

library(bsvars)


load("spartan/results/tax06nPM.rda")
post06 = post

load("spartan/results/taxMRnPM.rda")
postMR = post

load("spartan/results/tax23nPM.rda")
post23 = post

data("us_fiscal_lsuw")

rm("post","sddr")

# colors
cols    = c("darkorchid1",
            adjustcolor("darkorchid1", alpha.f = 0.3), 
            "orchid4",
            adjustcolor("orchid4", alpha.f = 0.3), 
            "darkorchid4",
            adjustcolor("darkorchid4", alpha.f = 0.3)
)

# the prior
#######################################################
set.seed(1)
n           = 10000
grid        = seq(from = -1.1, to = 1.1, by = 0.01)
s_          = spec$prior$sv_s_ / rchisq(n, df = spec$prior$sv_a_)
sigma2_omega = sapply(1:n, function(s){rgamma(1, shape = 1, scale = s_[s])})
omega_dprior = matrix(NA, n, length(grid))

for (i in 1:n) {
  omega_dprior[i,] = dnorm(x = grid, sd = sqrt(sigma2_omega[i]))
}
omega       = apply(omega_dprior, 2, mean)

# hist options
#######################################################
den_max = max(
  hist(post23$posterior$omega[1,], breaks = 200, plot = FALSE)$density,
  hist(post06$posterior$omega[1,], breaks = 200, plot = FALSE)$density,
  hist(postMR$posterior$omega[1,], breaks = 200, plot = FALSE)$density,
  hist(post23$posterior$omega[2,], breaks = 200, plot = FALSE)$density,
  hist(post06$posterior$omega[2,], breaks = 200, plot = FALSE)$density,
  hist(postMR$posterior$omega[2,], breaks = 200, plot = FALSE)$density,
  hist(post23$posterior$omega[3,], breaks = 200, plot = FALSE)$density,
  hist(post06$posterior$omega[3,], breaks = 200, plot = FALSE)$density,
  hist(postMR$posterior$omega[3,], breaks = 200, plot = FALSE)$density
)

main = c("2023 sample", "2006 sample", "MR sample")
xlab = c(expression(omega[ttr]), expression(omega[gs]), expression(omega[gdp]))
varr  = c("ttr", "gs", "gdp")

# histogram of omega_3(m)
#######################################################
pdf(
  file = "spartan/OMEGA_PM.pdf",
  height = 7,
  width = 12
)

par(
  mfcol = c(3,3),
  mar = c(4,5,1,0)
)
for (n in 1:3) {
  if (n > 1) {main = rep("", 3)}

  # 23 plot
  hist(
    post23$posterior$omega[n,], 
    col = cols[2], 
    border = cols[2], 
    breaks = 200,
    freq = FALSE,
    main = paste("volatility of volatility for", varr[n]),
    xlab = "",
    ylab = main[1],
    ylim = c(0, den_max),
    xlim = c(-1,1),
    cex.main = 1.6,
    cex.lab = 1.6,
    axes = FALSE
  )
  axis(2, 0:3, c(0,NA,NA,3), cex.axis = 1.6)
  axis(1, c(-1,0,1), c(-1,0,1), cex.axis = 1.6)
  lines(
    x = grid, 
    y = omega, 
    col = cols[1],
    lwd = 2
  )
  
  # 06 plot
  hist(
    post06$posterior$omega[n,], 
    col = cols[4], 
    border = cols[4], 
    breaks = 200,
    freq = FALSE,
    main = "",
    ylab = main[2],
    xlab = "",
    ylim = c(0, den_max),
    xlim = c(-1,1),
    cex.main = 1.6,
    cex.lab = 1.6,
    axes = FALSE
  )
  axis(2, 0:3, c(0,NA,NA,3), cex.axis = 1.6)
  axis(1, c(-1,0,1), c(-1,0,1), cex.axis = 1.6)
  lines(
    x = grid, 
    y = omega, 
    col = cols[3],
    lwd = 2
  )
  
  # MR plot
  hist(
    postMR$posterior$omega[n,], 
    col = cols[6], 
    border = cols[6], 
    breaks = 200,
    freq = FALSE,
    main = "",
    ylab = main[3],
    xlab = xlab[n],
    ylim = c(0, den_max),
    xlim = c(-1,1),
    cex.main = 1.6,
    cex.lab = 1.6,
    axes = FALSE
  )
  axis(2, 0:3, c(0,NA,NA,3), cex.axis = 1.6)
  axis(1, c(-1,0,1), c(-1,0,1), cex.axis = 1.6)
  lines(
    x = grid, 
    y = omega, 
    col = cols[5],
    lwd = 2
  )
}
dev.off()

system("cp spartan/OMEGA_PM.pdf /Users/twozniak/Research/Carnaval-gulasz-Fei/fei-tomasz-sv-prior/cg_paper/graph_OMEGA_PM.pdf")
