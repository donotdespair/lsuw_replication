
# DATA: data plot
############################################

data("us_fiscal_lsuw")
load("spartan/MR2006.RData")

data23    = us_fiscal_lsuw
data06    = ts(
              us_fiscal_lsuw[zoo::index(us_fiscal_lsuw) < 2007,], 
              start = c(1948, 1), frequency = 4
            )
dataMR    = ts(
              rbind(t(matrix(t(X)[1,1:12], 3,4)), t(Y)), 
              start = c(1950,1), frequency = 4
            )
dates23       = zoo::index(us_fiscal_lsuw)
dates06       = dates23[-((length(dates23) - 66):length(dates23))]
datesMR       = dates06[-(1:8)]
rm("us_fiscal_lsuw", "X", "Y", "IV")

date00        = 1980 # level standardisation date

# colors
cols    = c("darkorchid1",
            adjustcolor("darkorchid1", alpha.f = 0.3), 
            "orchid4",
            adjustcolor("orchid4", alpha.f = 0.3), 
            "darkorchid4",
            adjustcolor("darkorchid4", alpha.f = 0.3)
)





pdf(
  file = "spartan/DATA.pdf",
  height  = 3,
  width   = 12
)
par(
  mfrow = c(1,3),
  mar = c(4,4,1,1)
)
for (n in 1:3) {
  plot_range    = range(
    data23[,n] - data23[which(dates23 == date00),n],
    data06[,n] - data06[which(dates06 == date00),n],
    dataMR[,n] - dataMR[which(datesMR == date00),n]
  )
  plot(
    x = dates23,
    y = data23[,n] - data23[which(dates23 == date00),n],
    ylim = plot_range,
    type = "l",
    lwd = 2,
    col = cols[1],
    main = colnames(data23)[n],
    ylab = "",
    xlab = "time",
    bty = "n",
    axes = FALSE,
    cex.main = 1.6,
    cex.lab = 1.6
  )
  axis(1, cex.axis = 1.6,)
  axis(2, 
       c(plot_range[1], 0, plot_range[2]), 
       c(round(plot_range[1], 2), NA, round(plot_range[2], 2)), 
       cex.axis = 1.6,
  )
  lines(
    x = dates06,
    y = data06[,n] - data06[which(dates06 == date00),n],
    col = cols[3],
    lwd = 2
  )
  lines(
    x = datesMR,
    y = dataMR[,n] - dataMR[which(datesMR == date00),n],
    col = cols[5],
    lwd = 2
  )
}
dev.off()

system("cp spartan/DATA.pdf /Users/twozniak/Research/Carnaval-gulasz-Fei/fei-tomasz-sv-prior/cg_paper/graph_DATA.pdf")
