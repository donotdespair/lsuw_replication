
sddrs = matrix(NA, 0, 6)
for (iteration in 1:30) {
  load(paste0("JOE_RR_sim/code/results/bsvar_sv_",iteration,".rda"))
  sddrs_tmp = sddr
  load(paste0("JOE_RR_sim/code/results/bsvar_sv_5_",iteration,".rda"))
  sddrs_tmp = c(sddrs_tmp, sddr)
  load(paste0("JOE_RR_sim/code/results/bsvar_sv_chan_",iteration,".rda"))
  sddrs_tmp = c(sddrs_tmp, sddr)
  load(paste0("JOE_RR_sim/code/results/bsvar_sv_chan_5_",iteration,".rda"))
  sddrs_tmp = c(sddrs_tmp, sddr)
  load(paste0("JOE_RR_sim/code/results/bsvar_msh_",iteration,".rda"))
  sddrs_tmp = c(sddrs_tmp, sddr)
  load(paste0("JOE_RR_sim/code/results/bsvar_msh_5_",iteration,".rda"))
  sddrs_tmp = c(sddrs_tmp, sddr)
  sddrs = rbind(sddrs, sddrs_tmp)
}
colnames(sddrs) = c("sv", "sv_5", "chan", "chan_5", "msh", "msh_5")
apply(sddrs, 2, mean)
apply(sddrs, 2, sd)
apply(sddrs, 2, range)

# read spartan outputs
######################################################
B       = matrix(c(100, -20, 80, 200), 2, 2); B
phrase = c("sv_7", "sv_l", "sv_f", "sv_c")
process = c("hom","hom","het","het","hom","het","het","hom","het","het","hom","het","het")

files       = list.files("JOE_RR_sim/spartan/results/")
files       = files[grepl(phrase[4], files)]

sddrs       = rep(NA, length(files))
B_mean      = array(NA, c(2, 2, length(files)))
sigma_mean  = array(NA, c(2, 780, length(files)))
omega_mean  = matrix(NA, 2, length(files))

for (file in files) {
  load(paste0("JOE_RR_sim/spartan/results/",file))
  i               = which(file == files)
  sddrs[i]        = sddr
  B_mean[,,i]     = apply(post$posterior$B, 1:2, mean)
  omega_mean[,i]  = apply(abs(post$posterior$omega), 1, mean)
  sigma_mean[,,i] = apply(post$posterior$sigma, 1:2, mean)
}
names(sddrs) = rep(process, each = 2)
matrix(sddrs, ncol = 4)
apply(matrix(sddrs, ncol = 4), 2, \(x)(mean(x>0)))
B_mean
colnames(omega_mean) = rep(process, each = 2)
omega_mean
i = 1
plot.ts(t(sigma_mean[,,i])); i; i = i + 1

