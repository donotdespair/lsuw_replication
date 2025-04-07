
Rcpp::sourceCpp("src_cpp/bsvar_sv_chan.cpp")

options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
args
iteration <- as.integer(args[1])
iteration
rm(args)

library(bsvars)
set.seed(123 + iteration)
S_burn          = 1e4
S               = 1e4

T               = 260
B_dgp           = matrix(c(100, -20, 80, 200), 2, 2)
id              = paste0(paste0(rep(0, 4 - nchar(as.character(iteration))), collapse = ""), as.character(iteration))

yy              = as.matrix(read.csv(paste0("dgps/dgp", T, "_", id, ".csv"), header = FALSE))
N               = ncol(yy)
T               = nrow(yy)

B               = matrix(TRUE, N, N)
spec            = specify_bsvar_sv$new(yy, B = B) 
prior           = spec$prior$get_prior()
starting_values = spec$starting_values$get_starting_values()
starting_values$h = matrix(rnorm(prod(dim(yy)), sd = 0.01), ncol(yy), nrow(yy))
starting_values$S = matrix(1, ncol(yy), nrow(yy))
starting_values$B = B_dgp
starting_values$omega = c(0, sqrt(0.5))
starting_values$sigma2_omega = rep(10, N)
starting_values$rho = c(0, 0.9)
VB              = spec$identification$get_identification()

burn            = bsvar_sv_chan_cpp(S = S_burn, Y = t(yy), prior, VB, starting_values, thin = 1, centred_sv = FALSE, show_progress = TRUE)
post            = bsvar_sv_chan_cpp(S = S, Y = t(yy), prior, VB, starting_values = burn$last_draw, thin = 1, centred_sv = FALSE, show_progress = TRUE)

P     = pick_permutation_bsvars(post$last_draw$B, B_dgp, T, 
                                covariance_inv_bsvars(B_dgp, post$last_draw$sigma^2, T, communication_matrix(N, N)), 
                                allPermutations(N) - 1, 
                                allSigns(N), 
                                communication_matrix(N, N))
B_hat = P %*% post$last_draw$B
post$posterior = normalise_jaro_bsvar_sv (post$posterior, B_hat = B_hat)
# save(post, file = paste0("results/chan_post_", T, "_", id, ".rda"))
post$posterior$B[,,1]

sddr            = -as.numeric(verify_volatility_sv_chan_cpp( post$posterior, prior, t(yy) ))
save(sddr, file = paste0("results/s_chan_sddr_", T, "_", id, ".rda"))

load("dgps/dgp_B.rda")
Bm = mean((apply(post$posterior$B, 1:2, mean) - B)^2)

load(paste0("dgps/dgp_sds", T, "_", id, ".rda"))
sm = mean((apply(post$posterior$sigma, 1:2, mean) - t(sds))^2)

Sigma = array(0, c(N, N, T))
Sigma_tmp = matrix(NA, N, N)
for (s in 1:S) {
  Bi = solve(post$posterior$B[,,s])
  for (t in 1:T) {
    Sigma_tmp =  Bi %*% diag(post$posterior$sigma[,t,s]^2) %*% t(Bi)
    Sigma[,,t] =  Sigma[,,t] + Sigma_tmp
  }
}
Sigma = Sigma / S

SS = array(0, c(N, N, T))
BBi = solve(B)
for (t in 1:T) {
  SS[,,t] = BBi %*% diag(sds[t,]^2) %*% t(BBi)
}
Sm = mean((Sigma - SS)^2)

save(Bm, sm, Sm, file = paste0("results/s_chan_mse_", T, "_", id, ".rda"))
