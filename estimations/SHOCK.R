
# Check tax shocks identification
############################################################

# read the instrument
iv_tmp      = read.csv("data/MR2014-data/data-MR2014.csv")
iv          = xts::xts(iv_tmp$Tax.Narrative..TN., zoo::as.yearqtr(iv_tmp$Date))
iv          = iv[-(1:5)] 
iv_dates    = zoo::index(iv)

# read others' shocks
shock_bp    = read.csv("spartan/BP_reproductions/taxshocks_BP.csv")[,1]
shock_bp    = xts::xts(shock_bp, iv_dates)
shock_mr    = read.csv("spartan/BP_reproductions/taxshocks_MR.csv")[,1]
shock_mr    = xts::xts(shock_mr, iv_dates)
shock_le    = read.csv("spartan/BP_reproductions/taxshocks_LE.csv")[,1]
shock_le    = xts::xts(shock_le, iv_dates)

# read my shocks
load("spartan/results/tax06nPM.rda")
shock_lsuw_06   = bsvars::compute_structural_shocks(post)
shock_lsuw_06   = apply(shock_lsuw_06, 1:2, mean)
shock_lsuw_06   = shock_lsuw_06[,-(1:9)]
shock_lsuw_06   = xts::xts(t(shock_lsuw_06), iv_dates)

load("spartan/results/taxMRnPM.rda")
shock_lsuw_mr   = bsvars::compute_structural_shocks(post)
shock_lsuw_mr   = apply(shock_lsuw_mr, 1:2, mean)
shock_lsuw_mr   = shock_lsuw_mr[,-1]
shock_lsuw_mr   = xts::xts(t(shock_lsuw_mr), iv_dates)

load("spartan/results/tax23nPM.rda")
shock_lsuw_23   = bsvars::compute_structural_shocks(post)
shock_lsuw_23   = apply(shock_lsuw_23, 1:2, mean)
shock_lsuw_23   = shock_lsuw_23[,-(1:9)]
shock_lsuw_23   = shock_lsuw_23[,-(224:290)]
shock_lsuw_23   = xts::xts(t(shock_lsuw_23), iv_dates)

shock_corr = cbind(
  cor(cbind(iv, shock_lsuw_23))[2:4,1],
  cor(cbind(iv, shock_lsuw_06))[2:4,1],
  cor(cbind(iv, shock_lsuw_mr))[2:4,1],
  rbind(cor(cbind(iv, shock_bp, shock_mr, shock_le))[2:4,1], matrix(NA, 2, 3))
)
colnames(shock_corr) = c("2023-sample", "2006-sample", "MR-sample", "BP results", "MR results", "LE results")
rownames(shock_corr) = c("shock 1", "shock 2", "shock 3")

xtable::xtable(shock_corr, digits = 3)

