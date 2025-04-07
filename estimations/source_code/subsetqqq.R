sliceqq <-function(qqq,Subset=1){
  namelist = names(qqq)
  # qqqnew = qqq
  qqqnew = list()
  for (i in 1:length(namelist)) {
    if(length(dim(qqq[[namelist[i]]]))==3){
      dim = dim(qqq[[namelist[i]]])
      # qqqnew[[namelist[i]]]=qqq$posteriors[[namelist[i]]][,,Subset]
      qqqnew[[namelist[i]]] = array(qqq[[namelist[i]]][,,Subset],c(dim[1:2],length(Subset)))
    }else{
      qqqnew[[namelist[i]]]=qqq[[namelist[i]]][,Subset]
    }
    
  }
  return(qqqnew)
}