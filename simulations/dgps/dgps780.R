
#############################################
# SETUP
#############################################

# simulation setup
#############################################
S       = 100         # number of realisations for each dgp
T       = 780         # time dimension
N       = 2           # number of series
M       = 2           # number of regimes

# structural gdp parameters
#############################################
B       = matrix(c(100, -20, 80, 200), N, N); B # structural matrix
save(B, file = "JOE_RR_sim/spartan/dgps/dgp_B.rda")

Bi      = solve(B)
cov2cor(Bi %*% t(Bi))
# SV gdp parameters
#############################################
sv_rho  = 0.92         # persistence of the stochastic volatility
sv_s2   = 0.25        # stochastic volatility variance
sv_R    = diag(T)
mgcv::sdiag(sv_R, -1) = -sv_rho
sv_Ri   = solve(sv_R)

# GARCH gdp parameters
#############################################
ga_a    = 0.28
ga_b    = 0.7
ga_o    = 1 - ga_a - ga_b
ga_h0   = 1

# MSH gdp parameters
#############################################
ms_s2   = c(20, 10)
ms_PR   = matrix(c(.98, .02, .02, .98), M, M)
ms_pi   = c(1, 0)

#############################################
# HOMOSKEDASTIC
#############################################

# simulate data
# both shocks are homoskedastic
#############################################
initial_sim_no = 0
for (s in 1:S) {
  set.seed(1234 + s)
  
  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  sds     = matrix(1, T, N)            # standard deviations
  y       = (z * sds) %*% t(Bi)
  
  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y, 
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"), 
    sep = ",",
    row.names = FALSE, 
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

#############################################
# SV
#############################################

# simulate data
# first shocks is homoskedastic, second - SV
#############################################
initial_sim_no = S
for (s in 1:S) {
  set.seed(1234 + s)
  
  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  sds     = cbind(
    rep(1, T), 
    exp(0.5 * sv_Ri %*% matrix(rnorm(T), T, 1) * sqrt(sv_s2)) # conditional standard deviations
  )
  y       = (z * sds) %*% t(Bi)
  
  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y, 
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"), 
    sep = ",",
    row.names = FALSE, 
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}


# simulate data
# first shocks is SV, second - homoskedastic
#############################################
initial_sim_no = 2 * S
for (s in 1:S) {
  set.seed(1234 + s)
  
  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  sds     = cbind(
    exp(0.5 * sv_Ri %*% matrix(rnorm(T), T, 1) * sqrt(sv_s2)), # conditional standard deviations
    rep(1, T) 
  )
  y       = (z * sds) %*% t(Bi)
  
  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y, 
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"), 
    sep = ",",
    row.names = FALSE, 
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

# simulate data
# both shocks are SV
#############################################
initial_sim_no = 3 * S
for (s in 1:S) {
  set.seed(1234 + s)
  
  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  sds     = exp(0.5 * sv_Ri %*% matrix(rnorm(T * N), T, N) * sqrt(sv_s2)) # conditional standard deviations
  y       = (z * sds) %*% t(Bi)
  
  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y, 
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"), 
    sep = ",",
    row.names = FALSE, 
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

#############################################
# GARCH
#############################################

# simulate data
# first shocks is homoskedastic, second - GARCH
#############################################
initial_sim_no = 4 * S
for (s in 1:S) {
  set.seed(1234 + s)

  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  hs      = matrix(1, T, N)
  n       = 2
  for (t in 2:T) {
    hs[t, n] = ga_o + (ga_a * z[t - 1, n]^2 + ga_b) * hs[t - 1, n]
  }
  sds     = sqrt(hs)
  y       = (z * sds) %*% t(Bi)

  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y,
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}


# simulate data
# first shocks is GARCH, second - homoskedastic
#############################################
initial_sim_no = 5 * S
for (s in 1:S) {
  set.seed(1234 + s)

  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  hs      = matrix(1, T, N)
  n       = 1
  for (t in 2:T) {
    hs[t, n] = ga_o + (ga_a * z[t - 1, n]^2 + ga_b) * hs[t - 1, n]
  }
  sds     = sqrt(hs)
  y       = (z * sds) %*% t(Bi)

  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y,
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

# simulate data
# both shocks are GARCH
#############################################
initial_sim_no = 6 * S
for (s in 1:S) {
  set.seed(1234 + s)

  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  hs      = matrix(1, T, N)
  for (n in 1:N) {
    for (t in 2:T) {
      hs[t, n] = ga_o + (ga_a * z[t - 1, n]^2 + ga_b) * hs[t - 1, n]
    }
  }
  sds     = sqrt(hs)
  y       = (z * sds) %*% t(Bi)

  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y,
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

#############################################
# MSH
#############################################

# simulate data
# first shocks is homoskedastic, second - MSH
#############################################
initial_sim_no = 7 * S
for (s in 1:S) {
  set.seed(1234 + s)

  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  sds     = matrix(1, T, N)
  ms_homo = rep(1, T)
  n = 2
  count2 = 1
  while (count2 < 3) {
    for (t in 1:T) {
      if (t == 1) {
        ms_homo[t]  = sample(1:2, 1, prob = t(ms_PR) %*% ms_pi)
      } else {
        ms_homo[t]  = sample(1:2, 1, prob = t(ms_PR) %*% diag(M)[ms_homo[t - 1],])
      }
      if (ms_homo[t] == 2) {
        sds[t, n]   = sqrt(ms_s2[n])
      }
    }
    count2 = min(apply(as.matrix(sds[,n]) != 1, 2, sum))
  }
  y       = (z * sds) %*% t(Bi)

  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y,
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

# simulate data
# first shocks is MSH, second - homoskedastic
#############################################
initial_sim_no = 8 * S
for (s in 1:S) {
  set.seed(1234 + s)

  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  sds     = matrix(1, T, N)
  ms_homo = rep(1, T)
  n = 1
  count2 = 1
  while (count2 < 3) {
    for (t in 1:T) {
      if (t == 1) {
        ms_homo[t]  = sample(1:2, 1, prob = t(ms_PR) %*% ms_pi)
      } else {
        ms_homo[t]  = sample(1:2, 1, prob = t(ms_PR) %*% diag(M)[ms_homo[t - 1],])
      }
      if (ms_homo[t] == 2) {
        sds[t, n]   = sqrt(ms_s2[n])
      }
    }
    count2 = min(apply(as.matrix(sds[,n]) != 1, 2, sum))
  }
  y       = (z * sds) %*% t(Bi)

  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y,
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}


# simulate data
# both shocks are MSH
#############################################
initial_sim_no = 9 * S
for (s in 1:S) {
  set.seed(1234 + s)

  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  sds     = matrix(1, T, N)
  ms_homo = rep(1, T)
  n = 1:2
  count2 = 1
  while (count2 < 3) {
    for (t in 1:T) {
      if (t == 1) {
        ms_homo[t]  = sample(1:2, 1, prob = t(ms_PR) %*% ms_pi)
      } else {
        ms_homo[t]  = sample(1:2, 1, prob = t(ms_PR) %*% diag(M)[ms_homo[t - 1],])
      }
      if (ms_homo[t] == 2) {
        sds[t, n]   = sqrt(ms_s2[n])
      }
    }
    count2 = min(apply(as.matrix(sds[,n]) != 1, 2, sum))
  }
  y       = (z * sds) %*% t(Bi)

  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y,
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

#############################################
# HMSH
#############################################

# simulate data
# both shocks are HMSH
#############################################
initial_sim_no = 10 * S

for (s in 1:S) {
  set.seed(1234 + s)

  z       = matrix(rnorm(T * N), T, N) # standardised shocks
  n       = 1:2
  count2  = 1
  while (count2 < 3) {
    sds     = matrix(1, T, N)
    ms_homo = matrix(1, T, 2)
    for (i in n) {
      for (t in 1:T) {
        if (t == 1) {
          ms_homo[t,i]  = sample(1:2, 1, prob = t(ms_PR) %*% ms_pi)
        } else {
          ms_homo[t,i]  = sample(1:2, 1, prob = t(ms_PR) %*% diag(M)[ms_homo[t - 1,i],])
        }
        if (ms_homo[t,i] == 2) {
          sds[t, i] = sqrt(ms_s2[i])
        }
      }
    }
    count2 = min(apply(as.matrix(sds[,n]) != 1, 2, sum))
  }
  y       = (z * sds) %*% t(Bi)

  id_tmp  = s + initial_sim_no
  id      = paste0(paste0(rep(0, 4 - nchar(as.character(id_tmp))), collapse = ""), as.character(id_tmp))
  write.table(
    y,
    file = paste0("JOE_RR_sim/spartan/dgps/dgp", T, "_", id, ".csv"),
    sep = ",",
    row.names = FALSE,
    col.names = FALSE
  )
  save(sds, file = paste0("JOE_RR_sim/spartan/dgps/dgp_sds", T, "_", id, ".rda"))
}

