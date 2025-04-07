
# Estimation of the MS-SVAR-SV model with 3-level hierarchy for shrinkage
# 03: constant, and quadratic trend  in deterministic terms
# 06: as above but lower-triangular

# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

model = "06"

library(bsvars)
library(bsvarTVPs)

data("us_fiscal_lsuw")
data("us_fiscal_ex")
y           = us_fiscal_lsuw
T_tmp       = nrow(y)

p           = 4

Y           = t(y[(p + 1):nrow(y),])
X           = matrix(NA, nrow(y) - p, 0)
for (i in 1:p) {
  X         = cbind(X, y[(p + 1):nrow(y) - i,])
}
X           = t(cbind(X, 1, us_fiscal_ex[-(1:4),1:2]))
dates       = zoo::index(y)[(p + 1):nrow(y)]

T     = ncol(Y)
N     = nrow(Y)
K     = nrow(X)
M     = 2

B     = array(0, c(N, N, M))
for (m in 1:M) B[,,m] = m * diag(N)

set.seed(1234)
starting_values = list(
  B       = B,
  A       = cbind(1 * diag(N), matrix(0, N, K - N)),
  hyper   = matrix(1, 2 * N + 1, 2),
  PR_TR   = matrix(0.05, M, M) + 0.9*diag(M),
  xi      = rbind(
    c(rep(1, 150), rep(0, T - 150)),
    c(rep(0, 150), rep(1, T - 150))
  ),
  pi_0    = rep(1/M, M),
  h             = matrix(rnorm(N * T, sd = .01), N, T),
  rho           = rep(.5, N),
  omega         = matrix(.1, N, M),
  S             = matrix(1, N, T),
  sigma2_omega  = rep(1, N),
  s_            = rep(0.05, N)
)

VB     = list()
for (n in 1:N) {
  # VB[[n]]       = diag(N)
  VB[[n]]       = matrix(diag(N)[1:n,], ncol = N)
}

prior = list(
  B_V_inv    = diag(N),
  B_nu       = N,
  A          = cbind(diag(N), matrix(0, N, K - N)),
  A_V_inv    = diag(c(kronecker((1:p)^2, rep(1,N) ), rep(.01, K - p * N) )),
  hyper_nu_B = 10,
  hyper_a_B  = 10,
  hyper_s_BB = 100,
  hyper_nu_BB = 1,
  hyper_nu_A = 10,
  hyper_a_A  = 10,
  hyper_s_AA = 10,
  hyper_nu_AA = 10,
  PR_TR      = matrix(1, M, M) + 11*diag(M),
  sv_a_      = 1,
  sv_s_      = 0.1
)

set.seed(123)
system.time(
  qq0        <- bsvarTVPs::bsvar_mss_sv_boost(3e5, Y, X, prior, VB, starting_values, 1e3)
)

system.time(
  qq         <- bsvarTVPs::bsvar_mss_sv_boost(12e4, Y, X, prior, VB, qq0$last_draw, 2)
)

save(qq, prior, y, Y, X, VB, file = paste0("results/tax23ms.rda"))

qqq         = bsvarTVPs::normalise(qq, VB = VB)
bbb         = bsvarTVPs::structural_to_array(qqq)

save(qqq, bbb, prior, y, Y, X, VB, file = paste0("results/tax23ms", model,".rda"))
