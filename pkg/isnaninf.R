isnaninf=function(M){
  ini=any(is.nan(M)) | any(is.infinite(M))
  return(ini)
}