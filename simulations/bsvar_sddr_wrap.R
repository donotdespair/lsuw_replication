
# read spartan outputs
######################################################

T           = 780

S           = 100
phrases     = c(paste0(c("sv", "chan", "msh"),"_sddr_780"), paste0(c("s_sv", "s_chan", "s_msh"),"_sddr_260"))
process     = c("hom","sv_hom","sv_het","sv_het","ga_hom","ga_het","ga_het","ms_hom","ms_het","ms_het","hms_het")
sddrs       = matrix(NA, length(phrases), S * length(process))

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
    
    sddrs[iteration, i]      <- sddr
  }
  
}

# l-value test c: logSDDR = 0
lv = array(
  apply(
    sddrs, 
    1, 
    function(x){apply(matrix(x > 0, ncol = length(process)), 2, mean)}
  ),
  dim = c(length(process), length(phrases)),
  dimnames = list(process, phrases)
)
# q-value test c(a) & a = 0.05
qv = array(
  apply(
    sddrs, 
    1, 
    function(x){
      ca = quantile(matrix(x, ncol = length(process))[,1], probs = 0.95)
      apply(matrix(x > ca, ncol = length(process)), 2, mean)
    }
  ),
  dim = c(length(process), length(phrases)),
  dimnames = list(process, phrases)
)

lv[,c(4,1,2,5)]

squares = list()
for (j in 1:4) {
  i = c(1,4,2,5)[j]
  sq_tmp = matrix(lv[2:10,i], ncol = 3)
  squares[[j]] = rbind(
    rep(lv[1,i], 4),
    cbind(
      sq_tmp,
      c(sq_tmp[1:2,3], lv[11,i])
    )
  )
  # squares[[j]][3:4,] = 1 - squares[[j]][3:4,]
  sq_tmp = matrix(qv[2:10,i], ncol = 3)
  squares[[4 + j]] = rbind(
    rep(qv[1,i], 4),
    cbind(
      sq_tmp,
      c(sq_tmp[1:2,3], qv[11,i])
    )
  )
  # squares[[4 + j]][3:4,] = 1 - squares[[4 + j]][3:4,]
}
squares
out = rbind(
  cbind(
    rbind(
      squares[[1]],
      squares[[2]]
    ),
    rbind(
      squares[[3]],
      squares[[4]]
    )
  ),
  cbind(
    rbind(
      squares[[5]],
      squares[[6]]
    ),
    rbind(
      squares[[7]],
      squares[[8]]
    )
  )
)
knitr::kable(
  out, 
  format = "latex",
  digits = 2
)



### Lütkepohl Woźniak 2020 
#############################################
squares = list()
for (j in 1:2) {
  i = c(3,6)[j]
  sq_tmp = matrix(lv[2:10,i], ncol = 3)
  squares[[j]] = rbind(
    rep(lv[1,i], 4),
    cbind(
      sq_tmp,
      c(sq_tmp[1:2,3], lv[11,i])
    )
  )
  # squares[[j]][3:4,] = 1 - squares[[j]][3:4,]
  sq_tmp = matrix(qv[2:10,i], ncol = 3)
  squares[[2 + j]] = rbind(
    rep(qv[1,i], 4),
    cbind(
      sq_tmp,
      c(sq_tmp[1:2,3], qv[11,i])
    )
  )
  # squares[[4 + j]][3:4,] = 1 - squares[[4 + j]][3:4,]
}
squares
out = rbind(
  rbind(
    squares[[1]],
    squares[[2]]
  ),
  rbind(
    squares[[3]],
    squares[[4]]
  )
)
knitr::kable(
  out, 
  format = "latex",
  digits = 2
)

