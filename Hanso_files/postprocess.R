postprocess=function(x,g,dnorm,X,G,w){
  for(j in 1:length(X)) dist[j]=norm(x-X[,j])
  evaldist=max(dist)
  mindist=min(dist)
  indx=which.min(dist)
  if(mindist==0 & indx==1)
  else if (mindist==0 & index >1){
    X[,c(1,indx)]=X[,c(indx,1)]
    G[,c(1,indx)]=G[,c(indx,1)]
    w[c(1,indx)]=w[c(1,indx)]
  }
  else{
    X=cbind(x,X)
    G=cbind(g,G)
    tmp=qpspecial(G)
    w=tmp$x
    d=tmp$d
    dnorm=norm(d)
  }
  loc$dnorm=dnorm
  loc$evaldist=evaldist
  return(lsit(loc,X,G,w))
}