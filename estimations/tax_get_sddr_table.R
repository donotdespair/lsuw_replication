
# bsvars package installed from the developer's repo using:
# devtools::install_git("https://github.com/bsvars/bsvars.git")

files_in_results  = list.files("spartan/results/")
files_sddr        = files_in_results[grepl("sddr_", files_in_results)]
I                 = length(files_sddr)
sddrs             = array(
                      NA, 
                      c(3,2,3,12), 
                      dimnames = list(c("ttr", "gs", "gdp"), 
                                      c("logSDDR", "log_SDDR_se"), 
                                      c("23", "06", "MR"), 
                                      1:12)
                      )
models            = c("23", "06", "MR")

for (n in 1:length(models)) {
  files_sddr_tmp    = files_sddr[grepl(models[n], files_sddr)]
  
  for (i in 1:12) {
    load(paste0("spartan/results/", files_sddr_tmp[i]))
    sddrs[,, n, i]   = cbind(sddr$logSDDR, sddr$log_SDDR_se)
  }
}

apply(sddrs, 1:3, mean)
apply(sddrs, 1:3, sd)

range(sddrs[,1,1,]) # the same conclusions
range(sddrs[,1,2,]) # the same conclusions (except for gs shock for model 12)
range(sddrs[,1,3,]) # the same conclusions
