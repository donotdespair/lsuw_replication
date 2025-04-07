
# PAR: structural form parameter estimates
############################################

library(bsvars)

source("spartan/source_code/TransformtoMR.R")
source("spartan/source_code/cut_sf_MCMC.R")
source("spartan/source_code/bracket.R")

load("spartan/results/tax23.rda")
post23 = post

load("spartan/results/tax06n.rda")
post06 = post

load("spartan/results/taxMRn.rda")
postMR = post

data("us_fiscal_lsuw")

rm("post","sddr","spec")

sf23    = cut_sf_MCMC(TranstoMR(post23$posterior))
sf06    = cut_sf_MCMC(TranstoMR(post06$posterior))
sfMR    = cut_sf_MCMC(TranstoMR(postMR$posterior))

ord       = c(3,4,1,2,5,6,7,8,9)
par_table = cbind(
  round(apply(sf23, 2, mean)[ord], 3), bracket(apply(sf23, 2, sd)[ord]),
  round(apply(sf06, 2, mean)[ord], 3), bracket(apply(sf06, 2, sd)[ord]),
  round(apply(sfMR, 2, mean)[ord], 3), bracket(apply(sfMR, 2, sd)[ord])
)
xtable::xtable(par_table)
