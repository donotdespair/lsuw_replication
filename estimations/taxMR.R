
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

set.seed(1234)
B         = matrix(TRUE, 3, 3)
spec      = specify_bsvar_sv$new(
              data = yMR,
              p    = 4,
              exogenous = exMR,
              B    = B
            )

burn      = estimate(spec, 3e5, thin = 1e4)
post      = estimate(burn, 6e5, thin = 10)

sddr      = verify_volatility(post)

save(spec, post, sddr, file = paste0("results/taxMR.rda"))
