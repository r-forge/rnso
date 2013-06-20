mat1=matrix(c(0.5377,0.5381),1,2)
G=mat1
maxit=100
x=c()
m=dim(G)[1]
n=dim(G)[2]
#   if(m*n==0){
#     warning("G empty")
#     Return(list(x=c(),d=c(),q=Inf,info=c(2,0))
#   }
#maxit=max(varargin[[1]],10)  
e=rep(1,n)
#x=varargin[[2]]
nx=length(x)
#if(nx==0) {
#  x=as.matrix(e)
#}
x=e
#if(nx!=n) check error any(x)<0?
idx=seq(1,(n*n),by=n+1)
Q=t(G)%*%G
z=x
y=0
eta=0.9995
delta=3
mu0=as.numeric((t(x)%*%z)/n)
#print(mu0)
tolmu=1e-5
tolrs=1e-5
kmu=tolmu*mu0
print(kmu)
nQ=norm(Q,"I")+2
krs=tolrs*nQ
ap=0
ad=0
#for(k in 1:maxit){
  k=1
  r1=-Q%*%x+e*y+z
  r2=-1+sum(x)
  r3=-x*z
  rs=norm(rbind(r1,r2),"I") #check norm(,inf)?
  print(r3)
  (mu=-sum(r3)/n)
  
  if(as.numeric(mu)<as.numeric(kmu)){
    if(rs<krs) {info=c(0,k-1); break}
  }
  zdx=z/x
  QD=Q
  QD[idx]=QD[idx]+zdx
  C=chol(QD)
  KT=solve(t(C))%*%e
  M=t(KT)%*%KT
  r4=r1+r3/x
  r5=t(KT)%*%(solve(t(C))%*%r4)
  r6=r2+r5
  dy=-r6/M
  r7=r4+e*dy
  dx=solve(C)%*%(solve(t(C))%*%r7)
  dz=(r3-z*dx)/x
  p=-x/dx
  ap=min(min(p*(p>0)),1)
  ap=max(ap,1)
  p=-z/dz
  ad=min(min(p*(p>0)),1)
  ad=max(ad,1)
  mauff=(t(x+ap*dx)%*%(z+ad*dz))/n
  sig=(mauff/mu)^delta
  r3=r3+sig*mu
  r3=r3-dx*dz
  r4=r1+r3/x
  r5=t(KT)%*%(solve(t(C))%*%r4)
  r6=r2+r5
  dy=-r6/M
  r7=r4+e%*%dy
  dx=solve(C)%*%(solve(t(C))%*%r7)
  dz=(r3-z*dx)/x
  p=-x/dx
  ap=min(min(p*(p>0)),1)
  ap=max(ap,1)
  p=-z/dz
  ad=min(min(p*(p>0)),1)
  ad=max(ad,1)
  x=x+eta*ap*dx
  y=y+eta*ad*dy
  z=z+eta*ad*dz
#}
if(k==maxit) info=c(1,k)
x=max(x,0)
x=x/sum(x)
d=G*x
q=t(d)%*%d
retobj=list(x=x,d=d,q=q,info=info)