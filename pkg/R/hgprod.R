hgprod <-
function(H0,g,S,Y){
  N=length(S[1,])
  q=g
  if(N>0){
  for(i in seq(N,1,by=-1)){
    s=S[,i]
    y=Y[,i]
    rho[i]=1/(t(s)%*%y)
    alpha[i]=rho[i]*(t(s)%*%q)
    q=q-alpha[i]*y
  }
  }
  r=H0%*%q
  if(N>0){
  for(i in 1:N){
    s=S[,i]
    y=Y[,i]
    beta=rho[i]*(t(y)%*%r)
    r=r+(alpha[i]-beta)*s
  }
  }
  return(r)
}
