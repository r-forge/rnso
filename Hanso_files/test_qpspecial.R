mat1=matrix(c(0.5377,0.5381),1,2)
G=mat1
qpspecial(G)

trancone <- function(p){
  if(all(p<0)) return(1)
  else return(min(min(p[p>0]),1))
}