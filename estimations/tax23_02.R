# as tax23.R but with sv_s_ = 1

model = "_02"

# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

library(bsvars)
data("us_fiscal_lsuw")
data("us_fiscal_ex")

set.seed(1234)
B         = matrix(TRUE, 3, 3)
spec      = specify_bsvar_sv$new(
              data = us_fiscal_lsuw,
              p    = 4,
              exogenous = us_fiscal_ex,
              B    = B
            )

spec$prior$sv_s_ = 1

burn      = estimate(spec, 3e5, thin = 1e4)
post      = estimate(burn, 6e5, thin = 10)

sddr      = verify_volatility(post)

save(spec, post, sddr, file = paste0("results/tax23", model,".rda"))
