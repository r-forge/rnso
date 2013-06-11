gradsamp1run=function(x0,f0,g0,pars,options){
  samprad=options$samprad
  for(choice=1:length(samprad)){
    tmp=gradsampfixed(x0,f0,g0,samprad(choice),pars,options)
    if quitall return(list(x,f,g,dnorm,X,G,w))
    x0=x
    f0=f
    g0=g
  }
  return(list(x,f,g,dnorm,X,G,w))
}