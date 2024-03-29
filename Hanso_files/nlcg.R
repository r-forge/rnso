nlcg <- function(fn,gr,nvar,nstart=10,x0 = matrix(rnorm(nvar*nstart),nvar,nstart),
		     H0 = NULL, maxit = 1000,  fvalquit = -Inf,prtlevel=0,version="C",
                     normtol = 1e-6, xnormquit = Inf, evaldist = 1e-4, ngrad = 0,
                     scale = 1, wolfe1 = 1e-4, wolfe2 = 0.5, quitLSfail = TRUE,
		      strongwolfe = 1)
		      {
			x <- matrix(NA,nvar,nstart)
			f <- c()
			g <- matrix(NA,nvar,nstart)
			frec <- list()
			alpharec <- list()
			message <- list()
			for(run in 1:nstart){
			tmp <- nlcg1run(fn, gr, x0[,run], H0, maxit,  fvalquit,
                      strongwolfe,version,prtlevel,
                     normtol, xnormquit, evaldist, ngrad,
                     scale, wolfe1, wolfe2, quitLSfail)
			#print(tmp)
			x[,run] <- tmp$x
			f[run] <- tmp$f
			g[,run] <- tmp$g
			frec[[run]] <- tmp$frec
			alpharec[[run]] <- tmp$alpharec
			message[[run]] <- tmp$message
			}
			return(list(x = x, f = f, g = g, frec = frec, 
				alpharec = alpharec, message = message))
		      }