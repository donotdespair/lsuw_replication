# as tax06.R but with centred SV with Koop, Chan, Yu (2024) prior

# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

library(bsvars)
data("us_fiscal_lsuw")
data("us_fiscal_ex")

us_fiscal_ex    = ts(us_fiscal_ex[zoo::index(us_fiscal_ex) < 2007,], start = c(1948, 1), frequency = 4)
us_fiscal_lsuw  = ts(us_fiscal_lsuw[zoo::index(us_fiscal_lsuw) < 2007,], start = c(1948, 1), frequency = 4)

Rcpp::sourceCpp("src_cpp/bsvar_sv_ce.cpp")

set.seed(1234)
B         = matrix(TRUE, 3, 3)
spec      = specify_bsvar_sv$new(
  data = us_fiscal_lsuw,
  p    = 4,
  exogenous = us_fiscal_ex,
  B    = B,
  centred_sv = TRUE
)

Y               = spec$data_matrices$Y
X               = spec$data_matrices$X
prior           = spec$prior$get_prior()
VB              = spec$identification$get_identification()
starting_values = spec$starting_values$get_starting_values()

burn            = bsvar_sv_cpp(S = 3e5, Y, X, prior, VB, starting_values, thin = 1e4)
starting_values = burn$last_draw
post            = bsvar_sv_cpp(S = 6e5, Y, X, prior, VB, starting_values, thin = 10)
                               
save(spec, post, file = paste0("results/tax06_ce.rda"))
