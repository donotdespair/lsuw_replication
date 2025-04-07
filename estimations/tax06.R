
# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

library(bsvars)
data("us_fiscal_lsuw")
data("us_fiscal_ex")

us_fiscal_ex    = ts(us_fiscal_ex[zoo::index(us_fiscal_ex) < 2007,], start = c(1948, 1), frequency = 4)
us_fiscal_lsuw  = ts(us_fiscal_lsuw[zoo::index(us_fiscal_lsuw) < 2007,], start = c(1948, 1), frequency = 4)

set.seed(1234)
B         = matrix(TRUE, 3, 3)
spec      = specify_bsvar_sv$new(
              data = us_fiscal_lsuw,
              p    = 4,
              exogenous = us_fiscal_ex,
              B    = B
            )

burn      = estimate(spec, 3e5, thin = 1e4)
post      = estimate(burn, 6e5, thin = 10)

sddr      = verify_volatility(post)

save(spec, post, sddr, file = paste0("results/tax06.rda"))
