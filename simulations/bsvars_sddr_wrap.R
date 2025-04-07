
# read spartan outputs
######################################################

T           = 260

S           = 10
phrases     = paste0(c("s_sv", "s_chan", "s_msh"),"_sddr")
process     = c("hom","sv_hom","sv_het","sv_het","ga_hom","ga_het","ga_het","ms_hom","ms_het","ms_het","hms_het")
sddrs       = matrix(NA, length(phrases), S * length(process))

for (iteration in 1:3) {
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
    
    sddrs[iteration, i]      <- sddr
  }
  
}

array(
  apply(
    sddrs, 
    1, 
    function(x){apply(matrix(x > 0, ncol = length(process)), 2, mean)}
  ),
  dim = c(length(process), length(phrases)),
  dimnames = list(process, phrases)
)



# pos = matrix(sddrs[3,] > 0, ncol = length(process)); pos
