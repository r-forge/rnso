# old = 0
# new = 1
# fold = f0; gold = g0;
# xnew = x0 + new*d;
# [fnew,gradnew] = feval(fgname, xnew, pars);
# gnew = gradnew'*d;
# lszoom(old, new, ...
#           fold, fnew, gold, gnew, f0, g0, x0, d, pars, c1, c2, prtlevel);
# 
# x0 = c(-1.2,1)
# d=c(1,1)
# lo = 0
# hi = 1
# f0 = fr(x0)
# grad0 = grr(x0)
# flo = f0
# g0 <- sum(grad0*d)
# glo = g0
# xnew = x0+d
# fnew <- fr(xnew)
# gradnew <- grr(xnew)
# gnew <- sum(gradnew*d)
# flo2=fnew
# glo2=gnew
# cubic_interp(lo, hi, flo, flo2, glo, glo2)
# old <- 0
#   		fold <- f0
# 			gold <- g0
# 			new <- 1
#  xnew <- x0 + new*d
#   		  fnew <- fr(xnew)
# 			  gradnew <- grr(xnew)
# 			  gnew <- sum(gradnew*d)
# tmp <- lszoom(fr, grr, old, new, 
#   			    fold, fnew, gold, gnew, f0, g0,
# 				    x0, d, c1=0, c2=0.5, prtlevel=1)
# 
# 
# lo, hi, flo, fhi, glo, ghi, f0, g0, x0, d, pars, c1, c2, prtlevel

source('lszoom.R')
source('linesch_sw.R')
res=linesch_sw(fr,grr,c(-1.2,1),c(1,1))
