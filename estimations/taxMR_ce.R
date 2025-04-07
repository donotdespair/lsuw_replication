
# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

library(bsvars)
load("MR2006.RData")

yMR   = ts(
  rbind(
    t(matrix(t(X)[1,1:12], 3,4)),
    t(Y)
  ), start = c(1950,1), frequency = 4
)

tr = -(115.5:112.5)
exMR  = ts(
  rbind(
    cbind(tr, tr^2, rep(0,length(tr))),
    t(X)[,14:16]
  ), start = c(1950,1), frequency = 4
)

Rcpp::sourceCpp("src_cpp/bsvar_sv_ce.cpp")

set.seed(1234)
B         = matrix(TRUE, 3, 3)
spec      = specify_bsvar_sv$new(
  data = yMR,
  p    = 4,
  exogenous = exMR,
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
                               
save(spec, post, file = paste0("results/taxMR_ce.rda"))
