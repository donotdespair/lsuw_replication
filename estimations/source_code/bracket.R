bracket     = function(x, round_by = 3) {
  n         = length(x)
  pre       = rep("[", n)
  pos       = rep("]", n)
  paste0(paste0(pre,round(x,round_by)),pos)
}