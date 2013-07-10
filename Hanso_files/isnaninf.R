isnaninf=function(M){
  ini <- any(is.na(M)) || any(is.infinite(M))
  return(ini)
}