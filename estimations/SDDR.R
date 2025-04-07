# SDDR: a table
############################################

library(bsvars)

load("spartan/results/tax23.rda")
sddr23 = sddr

load("spartan/results/tax06n.rda")
sddr06 = sddr

load("spartan/results/taxMRn.rda")
sddrMR = sddr

data("us_fiscal_lsuw")

bracket     = function(x) {
  n         = length(x)
  pre       = rep("[", n)
  pos       = rep("]", n)
  paste0(paste0(pre,round(x,2)),pos)
  
}
sddr_table  = cbind(
                round(sddr23$logSDDR, 2), bracket(sddr23$log_SDDR_se),
                round(sddr06$logSDDR, 2), bracket(sddr06$log_SDDR_se),
                round(sddrMR$logSDDR, 2), bracket(sddrMR$log_SDDR_se)
              )

rownames(sddr_table) = colnames(us_fiscal_lsuw)

xtable::xtable(sddr_table)


# MR normalisation (as opposed to BP above)
############################################

load("spartan/results/tax23nMR.rda")
sddr23 = sddr

load("spartan/results/tax06nMR.rda")
sddr06 = sddr

load("spartan/results/taxMRnMR.rda")
sddrMR = sddr

data("us_fiscal_lsuw")

bracket     = function(x) {
  n         = length(x)
  pre       = rep("[", n)
  pos       = rep("]", n)
  paste0(paste0(pre,round(x,2)),pos)
  
}
sddr_table  = cbind(
  round(sddr23$logSDDR, 2), bracket(sddr23$log_SDDR_se),
  round(sddr06$logSDDR, 2), bracket(sddr06$log_SDDR_se),
  round(sddrMR$logSDDR, 2), bracket(sddrMR$log_SDDR_se)
)

rownames(sddr_table) = colnames(us_fiscal_lsuw)

xtable::xtable(sddr_table)

# PM normalisation (as opposed to BP above)
############################################

load("spartan/results/tax23nPM.rda")
sddr23 = sddr

load("spartan/results/tax06nPM.rda")
sddr06 = sddr

load("spartan/results/taxMRnPM.rda")
sddrMR = sddr

data("us_fiscal_lsuw")

bracket     = function(x) {
  n         = length(x)
  pre       = rep("[", n)
  pos       = rep("]", n)
  paste0(paste0(pre,round(x,2)),pos)
  
}
sddr_table  = cbind(
  round(sddr23$logSDDR, 2), bracket(sddr23$log_SDDR_se),
  round(sddr06$logSDDR, 2), bracket(sddr06$log_SDDR_se),
  round(sddrMR$logSDDR, 2), bracket(sddrMR$log_SDDR_se)
)

rownames(sddr_table) = colnames(us_fiscal_lsuw)

xtable::xtable(sddr_table)
