qpspecial=function(G,varargin$x=rep(1,n)){
  m=dim(G)[1]
  n=dim(G)[2]
  if(m*n==0) {print("G empty") Return}
  e=rep(1,n)
  x=varargin$x
  nx=length(x)
  #if(nx!=n) check error any(x)<0?
  idx=seq(1,(n*n),by=n+1)
  Q=t(G)%*%G
  z=x
  y=0
  eta=0.9995
  delta=3
  mu0=(t(x)%*%z)/n
  tolmu=1e-5
  tolrs=1e-5
  kmu=tolmu*mu0
  nQ=norm(Q,inf)+2
  krs=tolrs*nQ
  ap=0
  ad=0
  for(i in 1:maxit){
    r1=-Q*x+e*y+z
    r2=-1+sum(x)
    r3=-x*z
    rs=norm(c(r1,r2),inf) #check norm(,inf)?
    mu=-sum(r3)/n
    if(mu<kmu){
      if(rs<krs) {info=c(0,k-1) break}
    }
    zdx=z/x
    QD=Q
    QD[idx]=QD[idx]+zdx
    C=chol(QD)
    KT=solve(t(C))%*%e
    M=t(KT)%*%KT
    r4=r1+r3/x
    r5=t(KT)*(solve(t(C))*r4)
    r6=r2+r5
    dy=-r6/M
    r7=r4+e*dy
    dx=solve(C)*(solve(t(C))*r7)
    dz=(r3-z*dx)/x
    p=-x/dx
    ad=min(min(p(p>0)),1)
    if(ad==0) ad=1
    mauff=(t(x+ap*dx)*(z+ad*dz))/n
    sig=(mauff/mu)^delta
    r3=r3+sig*mu
    r3=r3-dx*dz
    r4=r1+r3/x
    r5=t(KT)*(solve(t(C))*r4)
    r6=r2+r5
    dy=-r6/M
    r7=r4+e*dy
    dx=solve(C)*(solve(t(C))*r7)
    dz=(r3-z*dx)/x
    p=-x/dx
    ad=min(min(p(p>0)),1)
    if(ad==0) ad=1
    x=x+eta*ap*dx
    y=y+eta*ad*dy
    z=z+eta*ad*dz
  }
  if(k==maxit) info=c(1,k)
  x=max(x,0)
  x=x/sum(x)
  d=G*x
  q=t(d)%*%d
  retobj=data.frame(x=x,d=d,q=q,info=info)
  return(retobj)
}























