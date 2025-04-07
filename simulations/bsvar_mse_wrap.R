
# read spartan outputs
######################################################

S           = 100
phrases     = c(paste0(c("sv", "chan", "ce"),"_mse_780"), paste0(c("s_sv", "s_chan", "s_ce"),"_mse_260"))
process     = c("hom","sv_hom","sv_het","sv_het","ga_hom","ga_het","ga_het","ms_hom","ms_het","ms_het","hms_het")
Bs          = matrix(NA, length(phrases), S * length(process))
ss          = matrix(NA, length(phrases), S * length(process))
Ss          = matrix(NA, length(phrases), S * length(process))

for (iteration in 1:length(phrases)) {
  phrase = phrases[iteration]
  cat(phrase, "\n")
  
  files     = list.files("JOE_RR_sim/spartan/results/")
  files     = files[grepl(phrase, files)]
  
  for (file in files) {
    # cat(file, "\n")
    tt = try(
      load(paste0("JOE_RR_sim/spartan/results/",file)) # change the path
    )
    i = as.numeric(substr(file, nchar(file) - 7, nchar(file) - 4))
    if (inherits(tt, "try-error")) next
    
    Bs[iteration, i]      <- Bm
    ss[iteration, i]      <- sm
    Ss[iteration, i]      <- Sm
  }
  
}
# Bs
# ss
# Ss


quant = 90

bbss = array(
  apply(
    Bs, 
    1, 
    function(x){apply(matrix(x, ncol = length(process)), 2, 
                      function(y){mean(y[order(y)][1:quant])})}
    ),
  dim = c(length(process), length(phrases)),
  dimnames = list(process, phrases)
)
bbss = sqrt(bbss)
Bmse = round(cbind(bbss[,1:3] / bbss[,1], bbss[,4:6] / bbss[,4]), 3)


ssss = array(
  apply(
    ss, 
    1, 
    function(x){apply(matrix(x, ncol = length(process)), 2, 
                      function(y){mean(y[order(y)][1:quant])})}
  ),
  dim = c(length(process), length(phrases)),
  dimnames = list(process, phrases)
)
ssss = sqrt(ssss)
smse = round(cbind(ssss[,1:3] / ssss[,1], ssss[,4:6] / ssss[,4]), 3)


SSss = array(
  apply(
    Ss, 
    1, 
    function(x){apply(matrix(x, ncol = length(process)), 2, 
                      function(y){mean(y[order(y)][1:quant])})}
  ),
  dim = c(length(process), length(phrases)),
  dimnames = list(process, phrases)
)
SSss = sqrt(SSss)
Smse = round(cbind(SSss[,1:3] / SSss[,1], SSss[,4:6] / SSss[,4]), 3)

out_mse = round(
  cbind(
    c(
      Bmse[1:4,2],
      Bmse[1:4,5]
    ),
    c(
      smse[1:4,2],
      smse[1:4,5]
    ),
    c(
      Smse[1:4,2],
      Smse[1:4,5]
    )
  ), 2
)

knitr::kable(
  out_mse, 
  format = "latex",
  digits = 2
)
