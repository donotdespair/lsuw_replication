
options(echo = TRUE)
args <- commandArgs(trailingOnly = TRUE)
args
iteration <- as.integer(args[1])
iteration
rm(args)

# read spartan outputs
######################################################

iteration   = 1
T           = 780

S           = 10
phrases     = c("sv", "chan", "msh")
process     = c("hom","sv_hom","sv_het","sv_het","ga_hom","ga_het","ga_het","ms_hom","ms_het","ms_het","hms_hom","hms_het","hms_het")
B           = matrix(c(100, -20, 80, 200), 2, 2)

sddrs       = array(NA, c(length(phrases), S * length(process)))
B_mean      = array(NA, c(2, 2, length(phrases), S * length(process)))
sigma_mean  = array(NA, c(2, T, length(phrases), S * length(process)))

# for (phrase in phrases[2:3]) {
phrase = phrases[iteration]
  cat(phrase, "\n")
  files     = list.files("results/")
  files     = files[grepl(phrase, files)]
  p         = which(phrase == phrases)
  
  for (file in files) {
    
    cat(file, "\n")
    tt = try(
      load(paste0("results/",file))
    )
    if (inherits(tt, "try-error")) next
    
    i                 = which(file == files)
    try(
      sddrs[p,i]      <- sddr
    )
    B_mean[,,p,i]     = apply(post$posterior$B, 1:2, mean)
    sigma_mean[,,p,i] = apply(post$posterior$sigma, 1:2, mean)
  }
# }
  save(sddrs, B_mean, sigma_mean, file = paste0("results/spartan_wrap",iteration,".rda")) 

# drs = apply(matrix(sddrs[3,] > 0, ncol = 13), 2, mean)
# names(drs) = process
# 
# B_rmse = rep(NA, length(process))
# sigma_rmse = rep(NA, length(process))
# names(B_rmse) = process
# names(sigma_rmse) = process
# for (i in 1:4) {
#   B_rmse[i]      = sqrt(mean((B_mean[,, 3, 1:10 + (i - 1) * 10] - array(rep(B, 10), c(2,2,10)))^2))
#   sigma_rmse[i]   = sqrt(mean((sigma_mean[,, 3, 1:10 + (i - 1) * 10] - 1)^2))
# }
# 
# save(sddrs, B_mean, sigma_mean, B_rmse, sigma_rmse, drs, file = paste0("results/spartan_wrap",iteration,".rda"))
