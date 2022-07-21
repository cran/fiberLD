## R routines for wood fiber length determination... 

# Fines and fibers Y densities...
fy.fines <- function(x,th.fines){
## gets density of fine length as density function of generalized gamma distribution (of library(VGAM)) 
## * x - vector of quantiles of true length of fine
## * th.fines - vector of genGamma parameters: (b.fines,d.fines,k.fines)
     ## * b.fines - scale parameter, b.fines >0 
     ## * d.fines, k.fines - positive parameters of gengamma distribution
  dgengamma.stacy(x=x,scale=th.fines[1],d=th.fines[2],k=th.fines[3])
}

fy.fibers <- function(x,th.fibers){
## gets density of fiber length as density function of generalized gamma distribution (of library(VGAM)) 
## * x - vector of quantiles of true length of fibers
## * th.fibers - vector of genGamma parameters: (b.fibers,d.fibers,k.fibers)
  dgengamma.stacy(x=x,scale=th.fibers[1],d=th.fibers[2],k=th.fibers[3])
}

# functions for get various length distributions...

p.uc <- function(y,r){
## gets probability for a cell of true length Y=y to be uncut in the increment core,
## p_uc(y)=P(X=y| Y=y) ...
## * y - the true cell (fine or fiber) length
## * r - radius of the increment core
 u = 2*r^2*asin((4*r^2-y^2)^.5/2/r)-y/2*(4*r^2-y^2)^.5
 t = pi*r^2+2*r*y
 z = u/t
 z[y>=2*r] = 0 
 z
}

fx.given.y <- function(x,y,r){
## gets density of X given Y=y (derivative of P(X<= y| Y=y) for x<y wrt x)...
## * x - observed cell (fine or fiber) length in the increment core
## * y - true cell length
## * r - radius of the core
  t = pi*r^2+2*r*y
  (8*r^2-3*x^2+x*y)/t/(4*r^2-x^2)^.5
}

# getting density of observed FINE length, fines X, ...

integrand.fines<-function(y,x,th.fines,r){
## integrand needed for calculating X fine density 
  fx.given.y(x,y,r)*fy.fines(y,th.fines)
}

fx.fines <- function(X,th.fines,r){
## gets density of observed fine length (X fines)... 
  z=c()
  for(i in 1:length(X)) {
      x=X[i]
      I=integrate(integrand.fines, lower=x, upper=100, x=x, th.fines=th.fines, 
         r=r, stop.on.error = F, rel.tol=.Machine$double.eps^.1)
      z[i]=p.uc(x,r)*fy.fines(x,th.fines)+I$value
  }
  z
} 

# getting density of observed FIBER length, fiber X, ...
integrand.fibers<-function(y,x,th.fibers,r){
## integrand needed for calculating X fiber density 
  fx.given.y(x,y,r)*fy.fibers(y,th.fibers)
}

fx.fibers <- function(X,th.fibers,r){
## gets density of observed fiber length (X fiber)... 
## * X - a vector of fiber length
   z=c()
   for(i in 1:length(X)){
      x=X[i]
      I=integrate(integrand.fibers,  lower=x, upper=100, x=x, th.fibers=th.fibers, 
           r=r, stop.on.error = F, rel.tol=.Machine$double.eps^.1)
      z[i]=p.uc(x,r)*fy.fibers(x,th.fibers)+I$value
   }
   z
} 

 fx.fines1 <- function(x,th.fines,r){
  ## gets density of a single observation of fine length...  
  ## * th.fines- vector of parameters (b.fines,d.fines,k.fines)
   I=integrate(integrand.fines, lower=x, upper=20, x=x, th.fines=th.fines, 
         r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
   z=p.uc(x,r)*fy.fines(x,th.fines=th.fines)+I$value
  } 

  fx.fibers1 <- function(x,th.fibers,r){
  ## gets density of a single observation of fiber length... 
    I=integrate(integrand.fibers, lower=x, upper=20, x=x, th.fibers=th.fibers,
           r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
    z=p.uc(x,r)*fy.fibers(x,th.fibers=th.fibers)+I$value
  } 

# Mixture density (with transformed parameters)...
fx.1 <-function(y,theta, r,log=TRUE){
## mixture density of a single value of the observed cell, fine and fiber...
## * theta - vector of seven transformed parameters(eps,b.fines,d.fines,k.fines,b.fibers,d.fibers,k.fibers) 
## * eps - a scalar, proportion of fines in the increment core
  theta[2:7] <- exp(theta[2:7])
  eps <- exp(theta[1])/(1+exp(theta[1])) 
 
  if(y>=(2*r) || y<=0) z <- 0
  else{
     z<- (1-eps)*fx.fibers1(y,th.fibers=theta[5:7],r) + eps*fx.fines1(y, th.fines=theta[2:4],r)
  }
  if (log) z <- log(z+1e-7) else z <- z+1e-7
  z
} ## fx.1


 lx <- function(x,theta,r, log=TRUE){
    ## vector of mixture densities of n cells length...
       sapply(x, FUN=fx.1, theta=theta,r=r,log=log) ## returns a vector
 }
  
 parlx <- function(cl,x,theta,r,log=TRUE){
      ## vector of mixture densities of n cells length...
       tmp <- parLapply(cl,x, fun=fx.1, theta=theta,r=r,log=log) 
       unlist(tmp)
 }  

## getting log likelihood with parallel computation....
loglik <- function(cl=NULL,theta,x, r){
## minus log likelihood function of generallized gamma model (mixture density)... 
## * theta - vector of transformed parameters     
   ## (eps,b.fines,d.fines,k.fines,b.fibers,d.fibers,k.fibers)
## * eps - a scalar, proportion of fines in the increment core
## * x - cell length
   if (is.null(cl))   
        ll <- sum(lx(x=x,theta=theta, r=r,log=TRUE))
   else       
      ll <- sum(parlx(cl=cl,x=x,theta=theta,r=r,log=TRUE))    
  -ll
} ## end loglik

 Dggamma.b <- function(y,th){
 ## gets derivative of gen gamma density wrt log(b)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      der <- d^2*b^(-d*k)*y^(d*k-1)*exp(-ybd)/gamma(k)*(ybd-k) ## wrt th_1
      der
 } 
 Dggamma.d <- function(y,th){
 ## gets derivative of gen gamma density wrt log(d)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      der <- d*b^(-d*k)*y^(d*k-1)*exp(-ybd)/gamma(k)*(1+log(ybd)*(k-ybd)) ## wrt th_2
      der
 } 
 Dggamma.k <- function(y,th){
 ## gets derivative of gen gamma density wrt log(k)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      der <- k*d*b^(-d*k)*y^(d*k-1)*exp(-ybd)/gamma(k)*(log(ybd)-digamma(k))  ## wrt th_3
      der
 } 


###############################################################################
## 2nd order partial derivatives of gen gamma density...
 D2ggamma.b2 <- function(y,th){
 ## gets direct 2nd derivative of gen gamma density wrt log(b)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      der <- d^3*b^(-d*k)*y^(d*k-1)*exp(-ybd)/gamma(k)*((ybd-k)^2-ybd ) 
      der
 } 
 D2ggamma.bd <- function(y,th){
 ## gets cross 2nd derivative of gen gamma density wrt log(b) log(d)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      dk <- d*k
      der <- d^2*b^(-dk)*y^(dk-1)*exp(-ybd)/gamma(k)*(-d*ybd*(ybd-k)*log(y/b)+
      d*ybd*log(y/b) + dk*(ybd-k)*log(y) - dk*log(b)*(ybd-k)+2*(ybd-k) ) 
      der
 } 
 D2ggamma.bk <- function(y,th){
 ## gets cross 2nd derivative of gen gamma density wrt log(b) log(k)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      dk <- d*k
      der <- dk*d*b^(-dk)*y^(dk-1)*exp(-ybd)/gamma(k)*(d*(ybd-k)*log(y)-
        digamma(k)*(ybd-k) - log(b)*d*(ybd-k) -1 ) 
      der
 } 
 D2ggamma.d2 <- function(y,th){
 ## gets direct 2nd derivative of gen gamma density wrt log(d)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
    #  aa <- (k-ybd)*log(y/b)
      ba <- 1+d*(k-ybd)*log(y/b) ## 1+d*aa
      der <- d*b^(-d*k)*y^(d*k-1)*exp(-ybd)/gamma(k)*( 2*ba-1 - ybd*log(ybd)^2 
         - ybd*log(ybd)*ba + k*log(ybd)*ba )
      der
 } 
 D2ggamma.dk <- function(y,th){
 ## gets cross 2nd derivative of gen gamma density wrt log(d) log(k)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      aa <- 1+d*(k-ybd)*log(y/b)
      der <- k* d*b^(-d*k)*y^(d*k-1)*exp(-ybd)/gamma(k)*( log(ybd)*aa -
               digamma(k)*aa  + log(ybd)  ) 
      der
 } 
 D2ggamma.k2 <- function(y,th){
 ## gets direct 2nd derivative of gen gamma density wrt log(k)...
      b <- th[1]; d <- th[2]; k <- th[3]
      ybd <- (y/b)^d
      dk <- d*k
      uu <- log(ybd)-digamma(k)
      der <- dk*b^(-dk)*y^(dk-1)*exp(-ybd)/gamma(k)*( k*log(ybd)*uu -
        k*digamma(k)*uu +  uu - k*trigamma(k)  )
      der
 } 

## checking 2nd  derivatives.....
#theta <- c(-.5,-2, -2,-3,-1,-2,-.5)

#theta <- c(0,log(.01),log(.1),log(10),rep(log(2),3)) 
#y <- 2.0
#th <- theta
#th[2:7] <- exp(th[2:7])
## eps <- exp(theta[1])/(1+exp(theta[1])) 
#d<- rep(NA,3)
#d[1] <-  D2ggamma.b2(y=y,th=th[2:4])
#d[2] <-  D2ggamma.bd(y=y,th=th[2:4])
#d[3] <-  D2ggamma.bk(y=y,th=th[2:4])
#del <- 1e-6
#fi <- rep(NA,3)
#for (i in 1:3){
# theta1 <- theta; theta1[i+1] <- theta[i+1]+del
# th1 <- theta1
# th1[2:7] <- exp(th1[2:7])
# fi[i] <- (Dggamma.b(y=y,th=th1[2:4])-Dggamma.b(y=y,th=th[2:4]))/del
#}
#(d-fi)/d


#d<- rep(NA,2)
#d[1] <-  D2ggamma.d2(y=y,th=th[2:4])
#d[2] <-  D2ggamma.dk(y=y,th=th[2:4])
#del <- 1e-6
#fi <- rep(NA,2)
#for (i in 1:2){
# theta1 <- theta; theta1[i+2] <- theta[i+2]+del
# th1 <- theta1
# th1[2:7] <- exp(th1[2:7])
# fi[i] <- (Dggamma.d(y=y,th=th1[2:4])-Dggamma.d(y=y,th=th[2:4]))/del
#}
#(d-fi)/d


#d <-  D2ggamma.k2(y=y,th=th[2:4])
#theta1 <- theta; theta1[4] <- theta[4]+del
#th1 <- theta1
#th1[2:7] <- exp(th1[2:7])
#fi <- (Dggamma.k(y=y,th=th1[2:4])-Dggamma.k(y=y,th=th[2:4]))/del
#(d-fi)/d



## end of 2nd order partial derivatives of gen gamma density
###################################################################


Dmixure.full <- function(x, theta,r){
## gets gradient of the LOG mixure distribution of the full model for a single observation wrt to all seven transformed para-s
## * theta - vector of transformed parameters
    ## (eps,b.fines,d.fines,k.fines,b.fibers,d.fibers,k.fibers)
    th <- theta
    theta[2:7] <- exp(theta[2:7])
    eps <- exp(theta[1])/(1+exp(theta[1])) 
    dmix <- rep(NA,7)
    dmix[1] <- eps^2/exp(th[1])*(fx.fines1(x=x,th.fines=theta[2:4],r=r) -
                  fx.fibers1(x=x,th.fibers=theta[5:7],r=r)) ## deriv wrt theta[1]

   
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*Dggamma.b(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
  
    dmix[2] <- eps*(p.uc(x,r=r)*Dggamma.b(y=x,th=theta[2:4]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*Dggamma.d(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
    dmix[3] <- eps*(p.uc(x,r=r)*Dggamma.d(y=x,th=theta[2:4]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*Dggamma.k(y,th), lower=x,upper=20, 
          x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
    dmix[4] <- eps*(p.uc(x,r=r)*Dggamma.k(y=x,th=theta[2:4]) + I$value)

    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*Dggamma.b(y,th), lower=x,upper=20, 
          x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
    dmix[5] <- (1-eps)*(p.uc(x,r=r)*Dggamma.b(y=x,th=theta[5:7]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*Dggamma.d(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
    dmix[6] <- (1-eps)*(p.uc(x,r=r)*Dggamma.d(y=x,th=theta[5:7]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*Dggamma.k(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1)
    dmix[7] <- (1-eps)*(p.uc(x,r=r)*Dggamma.k(y=x,th=theta[5:7]) + I$value)
    dmix/fx.1(y=x,theta=th, r=r,log=FALSE) 
} ## end Dmixure.full

 grad.mix.full <- function(x,theta,r){
  ## gradient matrix of mixture densities of n cells...
      sapply(x, FUN=Dmixure.full, theta=theta,r=r) ## returns a matrix
 } 

 pargrad.mix.full <- function(cl,x,theta,r){
  ## matrix of mixture densities of n cells length...
      tmp <- parLapply(cl,x, fun=Dmixure.full, theta=theta,r=r) 
      tmp <- unlist(tmp)
      matrix(tmp,length(theta),length(x))
 }  

dloglik <- function(cl,theta, x, r){
## gets the gradient of the log likelihood function of full mixure model...
## * theta - vector of transformed parameters 
   if (is.null(cl))    
       dll <- - rowSums(grad.mix.full(x=x,theta=theta,r=r))  
   else      
      dll <- -rowSums(pargrad.mix.full(cl=cl,x=x,theta=theta,r=r))     
   dll
} ## end dloglik


dloglik.nlm <- function(cl,theta, x, r){
## needed for nlm() function to get the loglik and its gradient as an attribute...
  ll <- loglik(cl,theta,x,r)
  attr(ll,"gradient") <- dloglik(cl,theta,x,r)
  ll
}



## checking derivatives by finite differences...
 ## checking deriv of log mixure for a single observ...
#theta <- c(-.5,-2, -2,-3,-1,-2,-.5)
#y <- 2.0
#dmix <- Dmixure.full(x=y,theta=theta,r=6)
#del <- 1e-6
#findif <- rep(0,7)
#for (i in 1:7){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (fx.1(y=y,theta=theta1, r=6,log=TRUE)-fx.1(y=y,theta=theta, log=TRUE))/del
#}
#(dmix-findif)/dmix
 ## checking deriv of log likelihood ...
#theta <- c(0,log(.01),log(.1),log(10),rep(log(2),3)) 
#y <- s[1:10]
#dll <- dloglik(cl=makeCluster(detectCores()),theta=theta,x=y,r=6)
#del <- 1e-5
#findif <- rep(0,7)
#for (i in 1:7){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik(cl=NULL,theta=theta1,x=y, r=6)-loglik(cl=NULL,theta=theta,x=y,r=6))/del
#}
#(dll-findif)/dll

D2mixure.full <- function(x, theta,r){
## gets Hessian of the LOG mixure distribution of the full model for a single observation wrt to all seven transformed para-s
## * theta - vector of transformed parameters
    ## (eps,b.fines,d.fines,k.fines,b.fibers,d.fibers,k.fibers)
    th <- theta
   # theta[2:7] <- exp(theta[2:7])
    theta <- exp(theta)
    eps <- theta[1]/(1+theta[1]) 
    h <- matrix(NA,7,7)
    grad <- Dmixure.full(x=x,theta=th,r=r) ## getting gradient of the log mixure
    h[1,1] <- -grad[1]^2 + (1- theta[1])*grad[1]/(1+theta[1])
    ind <- c(2:4)
    h[1,ind] <- h[ind,1] <- grad[ind]*(-grad[1] + eps/theta[1])
    ind <- c(5:7)
    h[1,ind] <- h[ind,1] <- -grad[1]*grad[ind] - eps^2/theta[1]*grad[ind]/(1-eps)
    efm <- eps/fx.1(y=x,theta=th, r=r,log=FALSE) 
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.b2(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[2,2] <- -grad[2]^2 + efm*(p.uc(x,r=r)*D2ggamma.b2(y=x,th=theta[2:4]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.bd(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[2,3] <- h[3,2] <- -grad[2]*grad[3] + efm*(p.uc(x,r=r)*D2ggamma.bd(y=x,th=theta[2:4]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.bk(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[2,4] <- h[4,2] <- -grad[2]*grad[4] + efm*(p.uc(x,r=r)*D2ggamma.bk(y=x,th=theta[2:4]) + I$value)
    ind <- c(5:7); h[2,ind] <- h[ind,2] <- -grad[2]*grad[ind]
    h[3,ind] <- h[ind,3] <- -grad[3]*grad[ind]
    h[4,ind] <- h[ind,4] <- -grad[4]*grad[ind]  
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.d2(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[3,3] <- -grad[3]^2 + efm*(p.uc(x,r=r)*D2ggamma.d2(y=x,th=theta[2:4]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.dk(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[3,4] <- h[4,3] <- -grad[3]*grad[4] + efm*(p.uc(x,r=r)*D2ggamma.dk(y=x,th=theta[2:4]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.k2(y,th), lower=x,upper=20, 
           x=x, th=theta[2:4],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[4,4] <- -grad[4]^2 + efm*(p.uc(x,r=r)*D2ggamma.k2(y=x,th=theta[2:4]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.b2(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    efm <- (1-eps)/fx.1(y=x,theta=th, r=r,log=FALSE) 
    h[5,5] <- -grad[5]^2 + efm*(p.uc(x,r=r)*D2ggamma.b2(y=x,th=theta[5:7]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.bd(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[5,6] <- h[6,5] <- -grad[5]*grad[6] + efm*(p.uc(x,r=r)*D2ggamma.bd(y=x,th=theta[5:7]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.bk(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[5,7] <- h[7,5] <- -grad[5]*grad[7] + efm*(p.uc(x,r=r)*D2ggamma.bk(y=x,th=theta[5:7]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.d2(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[6,6] <- -grad[6]^2 + efm*(p.uc(x,r=r)*D2ggamma.d2(y=x,th=theta[5:7]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.dk(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[6,7] <- h[7,6] <- -grad[6]*grad[7] + efm*(p.uc(x,r=r)*D2ggamma.dk(y=x,th=theta[5:7]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2ggamma.k2(y,th), lower=x,upper=20, 
           x=x, th=theta[5:7],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
    h[7,7] <- -grad[7]^2 + efm*(p.uc(x,r=r)*D2ggamma.k2(y=x,th=theta[5:7]) + I$value)
    h
} ## end D2mixure.full

## checking Hessian by finite differences...
 ## checking the Hessian of log mixure for a single observ...
#theta <- c(-.5,-2, -2,-3,-1,-2,-.5)
#y <- 2.0
#hh <- D2mixure.full(x=y,theta=theta,r=6)
#del <- 1e-5
#findif <- matrix(NA,7,7)
#for (i in 1:7) {
#   theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i,] <- (Dmixure.full(x=y,theta=theta1, r=6) - Dmixure.full(x=y,theta=theta,r=6 ))/del
#}
#(hh-findif)/hh


hess.mix.full <- function(x,theta,r){
  ## Hessian matrix of mixture densities of n cells...
      sapply(x, FUN=D2mixure.full, theta=theta,r=r) ## returns a matrix
 } 

 parhess.mix.full <- function(cl,x,theta,r){
  ## matrix of mixture densities of n cells length...
      tmp <- parSapply(cl,x, FUN=D2mixure.full, theta=theta,r=r) 
      tmp 
 }  

hessian <- function(cl,theta, x, r){
## gets the Hessian of the log likelihood function of full mixure model...
## * theta - vector of transformed parameters 
   if (is.null(cl))    
       dll <- matrix(rowSums(hess.mix.full(x=x,theta=theta,r=r)), 7,7)  
   else     
      dll <- matrix(rowSums(parhess.mix.full(cl=cl,x=x,theta=theta,r=r)),7,7)         
   dll
} ## end hessian



#################################################################
## getting log likelihood function for `Easy life',
## assuming uncut cells (generalized gamma distribution on X's)
## needed for parameters initialization....

fy.gg<-function(x,theta, log=TRUE){
## mixture density of cell length with genGamma on X's...
## * theta - vector of transformed parameters     
  theta[2:7] <- exp(theta[2:7])
  eps <- exp(theta[1])/(1+exp(theta[1])) 
  z=(1-eps)*fy.fibers(x = x,th.fibers = theta[5:7]) + eps*fy.fines(x = x,th.fines = theta[2:4])
  if (log) z <- log(z+1e-7) else z <- z+1e-7
  z
}
  

  fy <- function(x,theta, log=TRUE){
  ## mixture densities of n cells length...
     sapply(x, FUN=fy.gg, theta=theta,log=log) ## returns a vector
  }  

  parfy <- function(cl,x,theta, log=TRUE){
  ## mixture densities of n cells length...
     tmp <- parLapply(cl,x, fun=fy.gg, theta=theta,log=log) ## returns a vector
     unlist(tmp)
  }  

loglik.y.gg <- function(cl,theta,x){
## minus log likelihood function of generallized gamma model (mixture density) on X's... 
## * theta - vector of transformed parameters     
   ## (eps,b.fines,d.fines,k.fines,b.fibers,d.fibers,k.fibers)
## * eps - a scalar, proportion of fines in the increment core
## * x - cell length  
   if (is.null(cl))   
        ll <- sum(fy(x=x,theta=theta, log=TRUE))
   else       
      ll <- sum(parfy(cl=cl,x=x,theta=theta,log=TRUE))    
  -ll   
} ## end loglik.y.gg


 Dggamma <- function(y,th){
  ## gets gradient of the gen gamma density function wrt transformed parameters...
  ## * th - a vector of three untransformed parameters of gen gamma distribution
    b <- th[1]; d <- th[2]; k <- th[3]
    dggamma <- rep(NA,3)
    dk <- d*k
    ybd <- (y/b)^d
    p1 <- d*b^(-dk)*y^(dk-1)*exp(-ybd)/gamma(k)
    dggamma[1] <- p1*d*(ybd-k) ## wrt th_1 ## p1*d*(ybd-k)/b - derivative wrt b
    dggamma[2] <- p1*(1+log(ybd)*(k-ybd)) ## wrt th_2 ## p1*(1+log(ybd)*(k-ybd))/d - deriv wrt d
    dggamma[3] <- k*p1*(log(ybd)-digamma(k)) ## wrt th_3 ## p1*(log(ybd)-digamma(k)) - deriv wrt k
    dggamma 
  }



  Dmixure1 <- function(y,theta){
  ## gets gradient of the LOG mixure distribution for a single observation wrt to all seven transformed para-s
  ## * theta - vector of transformed parameters
    ## (eps,b.fines,d.fines,k.fines,b.fibers,d.fibers,k.fibers)
    th <- theta
    theta[2:7] <- exp(theta[2:7])
    eps <- exp(theta[1])/(1+exp(theta[1])) 
    dmix <- rep(NA,7)
    dmix[1] <- eps^2/exp(theta[1])*(fy.fines(x=y,th.fines=theta[2:4]) - fy.fibers(x=y,th.fibers = theta[5:7])) ## deriv w.r.t. theta1 (logit eps)
    dmix[2:4] <- eps*Dggamma(y=y, th=theta[2:4]) ## deriv w.r.t. b.fines, d.fines, k.fines
    dmix[5:7] <- (1-eps)*Dggamma(y=y, th=theta[5:7]) ## deriv w.r.t. b.fibers, d.fibers. k/fibers
   # dmix
    dmix/fy.gg(x=y,theta=th, log=FALSE)
  }

 grad.mix <- function(x,theta){
  ## matrix of gradient of mixture densities of n cells...
     sapply(x, FUN=Dmixure1, theta=theta) ## returns a matrix
  }  

  pargrad.mix <- function(cl,x,theta){
  ## matrix of gradient of mixture densities of n cells...
     tmp <- parLapply(cl,x, fun=Dmixure1, theta=theta) 
     tmp <- unlist(tmp)
     matrix(tmp,length(theta),length(x))  ## returns a matrix
  }  


loglik.y.grad <- function(cl,theta,x){
## calculates gradient of the log lik function of gen gamma mixure density on X's...
## * theta - vector of transformed parameters
  ## (eps,b.fines,d.fines,k.fines,b.fibers,d.fibers,k.fibers)  
  if (is.null(cl))    
       dll <- -rowSums(grad.mix(x=x,theta=theta))  
  else      
      dll <- -rowSums(pargrad.mix(cl=cl,x=x,theta=theta))
  dll  
}


## The end of `Easy life'
#################################

####################################################
## checking derivatives by finite differences ...
  ## checking deriv of gen gamma density...
#y <- 1.92
#th <- c(0.01,.1,10) 
#dg <- Dggamma(y=y,th=th) 
#theta <- log(th)
#del <- 1e-4
#findif <- rep(0,3)
#for (i in 1:3){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  th1 <- exp(theta1)
#  findif[i] <- (fy.fines(x=y,th.fines=th1)- fy.fines(x=y,th.fines=th))/del
#}
#(dg-findif)/dg  ## ok
  ## checking deriv of mixture wrt 7 para...
#theta <- c(-.5,-2, -2,-3,-1,-2,-.5)
#y <- 2.0
#dmix <- Dmixure1(y=y,theta=theta)
#del <- 1e-6
#findif <- rep(0,7)
#for (i in 1:7){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (fy.gg(x=y,theta=theta1, log=FALSE)-fy.gg(x=y,theta=theta, log=FALSE))/del
#}
#(dmix-findif)/dmix
## checking deriv of log mixure for a single observ...
#theta <- c(-.5,-2, -2,-3,-1,-2,-.5)
#y <- 2.0
#dmix <- Dmixure1(y=y,theta=theta)
#del <- 1e-6
#findif <- rep(0,7)
#for (i in 1:7){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (fy.gg(x=y,theta=theta1, log=TRUE)-fy.gg(x=y,theta=theta, log=TRUE))/del
#}
#(dmix-findif)/dmix

  ## checking deriv of log likelihood ...
#theta <- c(0,log(.01),log(.1),log(10),rep(log(2),3)) 
#y <- s[1:10]
#dll <- loglik.y.grad(theta=theta,x=y)
#del <- 1e-5
#findif <- rep(0,7)
#for (i in 1:7){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik.y.gg(theta=theta1,x=y)-loglik.y.gg(theta=theta,x=y))/del
#}
#(dll-findif)/dll


## end checking derivatives
##############################################

fled.kajaani <- function(data=stop("No data supplied"), r, model="ggamma", method="ML", parStart=NULL, fixed=NULL,optimizer=c("optim","L-BFGS-B","grad"), lower=-Inf, upper=Inf,cluster=0, ...) {
## function to estimate cell lengths from Kajaani data...
   
  if (!method%in%c("ML","SEM")) stop("unknown method for parameter estimation. 'ML' or 'SEM' methods are currently implemented")
 # if (!method%in%c("ML")) stop("unknown method for parameter estimation. Under the current version only 'ML' method is available")
  if (model=="ggamma" & method=="SEM")
            stop("'SEM' method is possible only with 'lognorm' model.")
  if (!(optimizer[1] %in% c("optim","nlm","nlm.fd")) )
          stop("unknown optimization method")
  if (model=="lognorm" & !is.null(fixed))
           stop("fixing parameters is possible only with 'ggamma' model.")
  if (optimizer[1]!="optim" & !is.null(fixed))
           stop("currently fixing parameters is available only with 'optim'.")
  userStart <- FALSE  ## no user-defined starting values of model parameters

  if (cluster==1){
    if (parallel::detectCores()>1)  ## no point otherwise
        cl <- parallel::makeCluster(parallel::detectCores()-1)
      else cl <- NULL
  } else if (cluster==0){
       cl <- NULL
  } else if (!is.null(cl)&&inherits(cluster,"cluster")){
     cl <- cluster
  } else {
       warning("Supplied cluster is unknown - ignored.")
       cl <- NULL
     } 

if (is.null(fixed) || sum(fixed)==0) { ## all parameters to be estimated
  if (!is.null(parStart)) {
      if (model=="ggamma") { 
           if (length(parStart)!=7) 
              stop("length of starting values of parameters, 'parStart', for ggamma model must be 7.")
           good <- parStart > 0
           if (sum(!good) > 0) stop("starting values of parameters must be positive")
           else if (parStart[1] > 1) stop("the starting value of the first parameter must be in (0,1).")
             else { userStart <- TRUE
                    thetaStart <- parStart
                    thetaStart[1] <- log(parStart[1]/(1-parStart[1]))
                    thetaStart[2:7] <- log(parStart[2:7])                 
              } 
      } else { ## if model =="lognorm"
               if (length(parStart)!=5) 
                  stop("length of starting values of parameters, 'parStart', for lognorm model must be 5.")
               if (parStart[3] < 0 | parStart[5] < 0)
                   stop("starting values of 3rd and 5th parameters must be positive") 
               else if (parStart[1] > 1 | parStart[1] < 0) stop("the starting value of the first parameter must be in (0,1).")
                  else {  userStart <- TRUE
                          thetaStart <- parStart
                          thetaStart[1] <- log(parStart[1]/(1-parStart[1]))
                          thetaStart[c(3,5)] <- log(parStart[c(3,5)])  
                   } 
        } ## if model =="lognorm"     
  } ## if (!is.null(parStart))
} else { ## if some para-s are fixed...
     if (length(fixed)!= 7) stop("length of 'fixed' must be 7.")
     else if (!is.logical(fixed)) stop("'fixed' must be of type 'logical'.")
     if (is.null(parStart)) stop("no fixed parameter values provided.")
     if (length(parStart)!=7) 
         stop("length of starting values of parameters, 'parStart' must be 7.")
     good <- parStart[fixed] > 0
     if (sum(!good)> 0) stop("fixed values of parameters must be positive")
     else if (fixed[1] && parStart[1]>1) stop("the fixed value of the first parameter must be in (0,1).")
     good <- parStart[!fixed] <= 0
     if (sum(!good) == sum(!fixed)) { 
                    userStart <- TRUE
                    thetaStart <- parStart
                    thetaStart[1] <- log(parStart[1]/(1-parStart[1]))
                    thetaStart[2:7] <- log(parStart[2:7]) 
     } else if (sum(good)< sum(!fixed)) 
            warning("only some of unfixed values of parameters are positive, user-specific starting values ignored.")  
} ## if some para-s are fixed
 
  if (model=="ggamma"){
       if (sum(is.finite(lower))==7){
           lower[1] <- log(lower[1]/(1-lower[1]))
           lower[2:7] <- log(lower[2:7])
       } else if (length(lower)>1 && length(lower)<7) stop("length of 'lower' must be 7.")
       if (sum(is.finite(upper))==7){
           upper[1] <- log(upper[1]/(1-upper[1]))
           upper[2:7] <- log(upper[2:7]) 
       } else if (length(upper)>1 && length(upper)<7) stop("length of 'upper' must be 7.")
  } else if (model=="lognorm"){
             if (sum(is.finite(lower))==5){
                  lower[1] <- log(lower[1]/(1-lower[1]))
                  lower[c(3,5)] <- log(lower[c(3,5)])
             } else if (length(lower)>1 && length(lower)<5) stop("length of 'lower' must be 5.")
             if (sum(is.finite(upper))==5){
                   upper[1] <- log(upper[1]/(1-upper[1]))
                   upper[c(3,5)] <- log(upper[c(3,5)])
             } else if (length(upper)>1 && length(upper)<5) stop("length of 'upper' must be 5.")
     } 

  x <- data

 # if (detectCores()>1)  ## no point otherwise
 #         cl <- makeCluster(detectCores()-1)
 # else cl <- NULL
  
  if (!userStart) {## get initial values of parameters ...
     if (!is.null(fixed) | sum(fixed)!=0){
          force(loglik.y.gg); force(fixed)
          thetaStart <- c(0,log(.01),log(.1),log(10),rep(log(2),3))
          .fixValues = parStart[fixed] 
          .thetaStart = thetaStart[!fixed]  ## initial values for para to be estimated 
          .loglik.y.gg <- function(cl,theta,...){ 
                 .theta = rep(NA,length(theta))  # rep(NA,sum(!fixed)) 
                 .theta[!fixed] = theta 
                 .theta[fixed] = .fixValues 
                 loglik.y.gg(cl,.theta,...) 
           }
          .loglik.y.grad <- function(cl,theta,...){ 
                 .gtheta = rep(NA,length(theta)) 
                 .gtheta[!fixed] = theta 
                 .gtheta[fixed] = .fixValues 
                 loglik.y.grad(cl,.gtheta,...)[!fixed]
           }
           .lower <- lower[!fixed]
           .upper <- upper[!fixed]
           b <- optim(par=.thetaStart,fn=.loglik.y.gg, gr=.loglik.y.grad, method="L-BFGS-B", x=x,cl=cl,lower=.lower, upper=.upper, ...)  
           fullpars <- rep(NA,7)
           fullpars[fixed] <- .fixValues 
           fullpars[!fixed]<- b$par
           b$par <- fullpars
     } else { if (model=="lognorm") {
                 loglik.y.gg <- loglik.y.logN
                 loglik.y.grad <- loglik.y.grad.logN
                 thetaStart <- c(0,-1,log(1.08),.6,log(.3))  ## rough start
               } else thetaStart <- c(0,log(.01),log(.1),log(10),rep(log(2),3)) ## rough start for gen gamma
              b <- optim(par=thetaStart,fn=loglik.y.gg, gr=loglik.y.grad, method="L-BFGS-B", x=x,cl=cl,lower=lower, upper=upper,...) ##lower=list(-2,-100,-3,-1,-1,-1,-1),upper=list(0,10,1,1,1,1,1)) 
            }
     thetaStart <- b$par
  } ## if (!userStart)

  
  par <- c() ## initializing estimated parameters
  object <- list()

if (method=="ML"){
  if (model=="lognorm") {
         loglik <- loglik.logN
         dloglik.nlm <- dloglik.logN.nlm
  }

  if (optimizer[1]=="optim"){
           if (!(optimizer[2] %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))){
                 warning("unknown optim() method, `Nelder-Mead' were used")
                  optimizer[2] <- "Nelder-Mead"
           }
           if (optimizer[2] %in% c("BFGS", "CG", "L-BFGS-B")){
               if (is.na(optimizer[3])){ 
                  warning("the third parameter of optimizer argument is not supplied, 
                     finite-difference approximation of the gradient were used")
                  grr <- NULL
               } else if (!(optimizer[3] %in% c("fd","grad"))){
                     warning("in the third parameter of optimizer, only `fd' and `grad'
                          options are possible, finite-difference 
                          approximation of the gradient were used")
                     grr <- NULL
               }  else if (optimizer[3] == "grad") { 
                               if (model=="lognorm") grr <- dloglik.logN
                               else grr <- dloglik
                       } else    grr <- NULL
           } else grr <- NULL

           ##      b <- optim(par=theta,fn=loglik, gr=grr, method=optimizer[2],x=x,r=r,cl=cl,...) 
           ## fixing ggamma para-s at the user-supplied values if needed ...
           if (is.null(fixed)) 
               b <- optim(par=thetaStart,fn=loglik, gr=grr, method=optimizer[2],x=x,r=r,cl=cl,lower=lower, upper=upper,...) 
           else{
             force(loglik)
             force(fixed)
             .fixValues = parStart[fixed] 
             .thetaStart = thetaStart[!fixed]  ## initial values for para to be estimated 
             .loglik <- function(cl,theta,...){ 
                 .theta = rep(NA,length(theta))  # rep(NA,sum(!fixed)) 
                 .theta[!fixed] = theta 
                 .theta[fixed] = .fixValues 
                 loglik(cl,.theta,...) 
              }
              if(!is.null(grr)){ 
                   .grr <- function(cl,theta,...){ 
                       .gtheta = rep(NA,length(theta)) 
                       .gtheta[!fixed] = theta 
                       .gtheta[fixed] = .fixValues 
                       grr(cl,.gtheta,...)[!fixed]
                    }
               } else{ .grr <- NULL }  
              .lower <- lower[!fixed]
              .upper <- upper[!fixed]
              b <- optim(par=.thetaStart,fn=.loglik, gr=.grr, method=optimizer[2],x=x,r=r,cl=cl,lower=.lower, upper=.upper, ...)  
              # object$fixed <- fixed
               object$estimated.par <- b$par
               par <- rep(NA,length(thetaStart))  ## full vector of paras
               par[fixed] <-.fixValues 
               par[!fixed] <- exp(b$par)
               if (!fixed[1]) par[1] <- exp(b$par[1])/(1+exp(b$par[1]))
          } ## !is.null(fixed)
     
           theta <- b$par
           if (model=="ggamma" && is.null(fixed)) {
                 par[1] <- exp(theta[1])/(1+exp(theta[1])) 
                 par[2:length(theta)] <- exp(theta[2:length(theta)])               
           } else if (model=="lognorm"){
                  par[1] <- exp(theta[1])/(1+exp(theta[1])) 
                  par[c(3,5)] <- exp(theta[c(3,5)]); par[c(2,4)] <- theta[c(2,4)] 
             } 
           
           ll <- b$value
           ## dll <- NULL
           iterations <- b$counts
           termcode <- b$convergence
           if (termcode == 0)
                     conv <- "Successful completion"
           else if (termcode == 1)  
                      conv <- "The iteration limit `maxit' had been reached"
           else if (termcode == 10)  
                      conv <- "Degeneracy of the Nelder-Mead simplex"
           else if (termcode == 51)  
                      conv <- "A warning from the `L-BFGS-B' method; see help for `optim' for further details"
           else if (termcode == 52)  
                      conv <- "An error from the `L-BFGS-B' method; see help for `optim' for further details"
      }   ## if (optimizer[1]=="optim")
       else if (optimizer[1]=="nlm"){ ## nlm() with analytical derivatives...
            b <- nlm(f=dloglik.nlm, p=thetaStart,iterlim=100,x=x,r=r,cl=cl) 
           } else if (optimizer[1]=="nlm.fd"){ ## nlm() with finite difference deriv-s...
                 b <- nlm(f=loglik, p=thetaStart,iterlim=100,x=x,r=r,cl=cl) 
             }

   if (optimizer[1]== "nlm.fd" || optimizer[1]== "nlm") {
         theta <- b$estimate
       #  par[2:7] <- exp(theta[2:7])
         par[1] <- exp(theta[1])/(1+exp(theta[1])) 
         if (model=="ggamma")
                 par[2:length(theta)] <- exp(theta[2:length(theta)])
         else {par[c(3,5)] <- exp(theta[c(3,5)]) ; par[c(2,4)] <- theta[c(2,4)] }
         ll <- b$minimum
         ## dll <- b$gradient
         iterations <- b$iterations
         termcode <- b$code
         if (termcode == 1)
                conv <- "Relative gradient is close to zero, current iterate is probably solution"
         else if (termcode == 2)  
                conv <- "Successive iterates within tolerance, current iterate is probably solution"
         else if (termcode == 3)  
                conv <- "Last global step failed to locate a point lower than `estimate'. Either `estimate' is an approximate local minimum of the function or `steptol' is too small"
         else if (termcode == 4)  
                conv <- "Iteration limit exceeded"
         else if (termcode == 5)  
                conv <- "Maximum step size `stepmax' exceeded five consecutive 
                   times. Either the function is unbounded below, becomes asymptotic 
                   to a finite value from above in some direction or stepmax is too small"   
    } ## if (optimizer[1]== "nlm.fd" || optimizer[1]== "nlm")

  object$loglik <- ll
  object$par <- par  ## estimated parameters on original scale
  if (model=="lognorm") names(object$par) <- c("eps","mu_fines", "sig_fines","mu_fibers", "sig_fibers")
  else names(object$par) <- c("eps","b_fines", "d_fines","k_fines","b_fibers", "d_fibers", "k_fibers")
  
  if (is.null(fixed))
     object$logpar <- theta ## estimated transformed parameters
  else {
     object$logpar <- par
     object$logpar[1] <- log(par[1]/(1-par[1]))
     object$logpar[2:7] <- log(par[2:7])
   } 

 
  ## object$dloglik <- dll 
  object$termcode <- termcode
  object$conv <- conv
  object$iterations <- iterations 
} ## end method=="ML"
   else if (method=="SEM"){
      kk.x <- function(x,r){
        kk <- (8*r^2-3*x^2+30*x)/((pi*r^2+x*r)*sqrt(4*r^2-x^2))
        kk
      }
      n <- length(x)
      n.iter <- 200; ## 200
      n.med <- 50;
      antal <- 100 ## 100 for mcmc
  
      if (!is.null(cl)) registerDoParallel(cl)
      rkoll <- p.uc(y=x,r=r)
      const <- kk.x(x=x,r=r)
      parStart <- thetaStart
      parStart[1] <- exp(thetaStart[1])/(1+exp(thetaStart[1])) 
      parStart[c(3,5)] <- exp(thetaStart[c(3,5)]); 
      parStart[c(2,4)] <- thetaStart[c(2,4)]
      th.fines <- parStart[2:3]
      th.fibers <- parStart[4:5]
      eps <- parStart[1]

      fin = x[x<0.5]
      fib = x[x>=0.5]
     n1 = length(fin);
      eps = n1/n
      th.fines <- c(mean(log(fin)),var(log(fin))^.5)
      th.fibers <- c(mean(log(fib)),var(log(fib))^.5)

      mu_fin_temp <- rep(0,n.iter)
      mu_fib_temp <- rep(0,n.iter)
      sigma_fin_temp <- rep(0,n.iter)
      sigma_fib_temp <- rep(0,n.iter)
      eps_temp <- rep(0,n.iter)
      randseeds1 <- matrix(sample(0:1e+5, n*n.iter,replace = TRUE), n.iter,n)
      randseeds2 <- matrix(sample(0:1e+5, n*n.iter,replace = TRUE), n.iter,n)

      for (i in 1:n.iter){ ## main EM loop...
        zvec <- rep(0,n)
        sumz <- 0
        sum_w.kv_fin <- matrix(0,n,2)
        sum_w.kv_fib <- matrix(0,n,2)

        fx_fin <- fyj.logN(x=x,th.j=th.fines)   
        fx_fib <- fyj.logN(x=x,th.j=th.fibers) 

        eps_diff <-1;
        while (eps_diff>1e-6){
           zvec <- eps*fx_fin/(eps*fx_fin+(1-eps)*fx_fib)
           sumz <- sum(zvec)
           eps_tem <- mean(zvec);
           eps_diff <- abs(eps-eps_tem);
           eps <- eps_tem;
        }   
     if (!is.null(cl)) {
          tmp <- foreach(j=1:n,  .combine='rbind') %dopar% {
                   tmp1 <- .Call("Simu1Cpp",x_=x[j], const_=const[j], mu_=th.fines[1], sigma_=th.fines[2], r_=r, r_koll_=rkoll[j], fx_fin_=fx_fin[j], zvec_=zvec[j], antal_=antal, rnd_seed_=randseeds1[i,j],PACKAGE="fiberLD") 
                   tmp2 <- .Call("Simu2Cpp",x_=x[j],const_=const[j], mu_=th.fibers[1], sigma_=th.fibers[2], r_=r, r_koll_=rkoll[j], fx_fib_=fx_fib[j], antal_=antal, rnd_seed_=randseeds2[i,j], PACKAGE="fiberLD")
                  c(tmp1, tmp2)
                  } 
          sum_w.kv_fin <- tmp[,1:2]
          sum_w.kv_fib <- tmp[,3:4] 
     } else{
             for (j in 1:n) {
             ## Fines Simulations
              sum_w.kv_fin[j,] <- .Call("Simu1Cpp",x_=x[j], const_=const[j], mu_=th.fines[1], sigma_=th.fines[2], r_=r, r_koll_=rkoll[j], fx_fin_=fx_fin[j], zvec_=zvec[j], antal_=antal, rnd_seed_=randseeds1[i,j],PACKAGE="fiberLD") 
             ## Fibers Simulation
             sum_w.kv_fib[j,] <- .Call("Simu2Cpp",x_=x[j],const_=const[j], mu_=th.fibers[1], sigma_=th.fibers[2], r_=r, r_koll_=rkoll[j], fx_fib_=fx_fib[j], antal_=antal, rnd_seed_=randseeds2[i,j], PACKAGE="fiberLD") 
            }
       }
  
         E1_fin <- sum(sum_w.kv_fin[,1]*zvec);
         E2_fin <- sum(sum_w.kv_fin[,2]*zvec);
         E1_fib <- sum(sum_w.kv_fib[,1]*(1-zvec));
         E2_fib <- sum(sum_w.kv_fib[,2]*(1-zvec));
    
         ## update parameters...
         mu_fin <- E1_fin/antal/sumz
         mu_fib <- E1_fib/antal/(n-sumz)
         sigma_fin <- (E2_fin/antal/sumz - mu_fin^2)^.5
         sigma_fib <- (E2_fib/antal/(n-sumz) - mu_fib^2)^.5
   
         ## save history...
         mu_fin_temp[i] <- mu_fin
         mu_fib_temp[i] <- mu_fib
         sigma_fin_temp[i] <- sigma_fin
         sigma_fib_temp[i] <- sigma_fib
         eps_temp[i] <- eps
      } ## end main EM loop
      mu_fin=mean(mu_fin_temp[(n.iter-n.med+1):n.iter])
      mu_fib=mean(mu_fib_temp[(n.iter-n.med+1):n.iter])
      sigma_fin=mean(sigma_fin_temp[(n.iter-n.med+1):n.iter])
      sigma_fib=mean(sigma_fib_temp[(n.iter-n.med+1):n.iter])
      ep=mean(eps_temp[(n.iter-n.med+1):n.iter])
      object$par <- par<- c(ep,mu_fin,sigma_fin, mu_fib, sigma_fib)
      names(object$par) <- c("eps","mu_fines", "sig_fines","mu_fibers", "sig_fibers")
      object$logpar <- object$par
      object$logpar[1] <- log(ep/(1-ep))
      object$logpar[c(3,5)] <- log(object$par[c(3,5)])
      object$loglik <- loglik.logN(cl=cl,theta=object$logpar,x=x, r=r)
   } ## end method=="SEM"

  ## getting things after fit ...
   if (model=="lognorm"){
      object$cov.logpar <- solve(-hessian.logN(cl=cl,theta=object$logpar, x=x, r=r))
      grad <- c(par[1]-par[1]^2, 1,par[3],1,par[5])
      object$cov.par <- diag(grad)%*%object$cov.logpar%*%diag(grad)
      # grad <- grad.mix.full.logN(x=x,theta=object$logpar, r=r)
      # object$outer.cov.par <- solve(grad%*%t(grad))
  } else{
         h <- hessian(cl=cl,theta=object$logpar, x=x, r=r)
         object$cov.logpar <- solve(-h)
         grad <- c(par[1]-par[1]^2, par[2:7])
         object$cov.par <- diag(grad)%*%object$cov.logpar%*%diag(grad)
      object$eigen <- eigen(-h)$values   
         if (!is.null(fixed)){ ## columns and rows for fixed parameters to be zeroed...
           ind.fix <- which(fixed) ## which indices are TRUE
           h <- h[-ind.fix,]
           h <- h[,-ind.fix]
           object$cov.logpar <- solve(-h)
       object$eigen <- eigen(-h)$values
        #  object$cov.logpar <- object$cov.logpar[-ind.fix,]
        #   object$cov.logpar <- object$cov.logpar[,-ind.fix]
           grad <- grad[-ind.fix] 
           object$cov.par <- diag(grad)%*%object$cov.logpar%*%diag(grad)    
         }
      # grad <- grad.mix.full(x=x,theta=object$logpar, r=r)
      # object$outer.cov.par <- solve(grad%*%t(grad))
    }

  if (!is.null(cl)) stopCluster(cl)
  ## getting estimates of the mean lengths of fiber and fine in a standing tree...
  if (model=="lognorm") {
      fy.fines <- fyj.logN
      th.fines <- par[2:3]
  } else  th.fines <- par[2:4]
  pir <- pi*r
  I=integrate(function(y,th.fines,pir) fy.fines(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
  object$mu.fines <- (1/I$value -pir)/2
  if (model=="lognorm") {
      fy.fibers <- fyj.logN
      th.fibers <- par[4:5]
  } else  th.fibers <- par[5:7]
  I=integrate(function(y,th.fibers,pir) fy.fibers(y,th.fibers)/(pir+2*y), lower=0,upper=20, th.fibers=th.fibers,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
  object$mu.fibers <- (1/I$value -pir)/2
  ## getting estimate of the expected value of the cell length in a standing tree, E(W)...
  object$mu.cell <- (2*object$mu.fines*object$mu.fibers+par[1]*pir*object$mu.fines+(1-par[1])*pir*object$mu.fibers)/(2*(par[1]*object$mu.fibers+(1-par[1])*object$mu.fines)+pir)   
  ## getting the proportion of fines in a standing tree...
  object$prop.fines <- par[1]*(pir+2*object$mu.cell)/(pir+2*object$mu.fines)
  object
} ## end fled.kajaani


#######################

fled.micro <- function(data=stop("No data supplied"), r, model="ggamma", parStart=NULL, optimizer=c("optim","L-BFGS-B","grad"), lower=-Inf, upper=Inf, ...) {
## function to estimate cell lengths from microscopy data...
  if (!(optimizer[1] %in% c("optim","nlm","nlm.fd")) )
          stop("unknown optimization method")
  
  userStart <- FALSE  ## no user-defined starting values of model parameters
  if (!is.null(parStart)) {
      if (model=="ggamma") { 
           if (length(parStart)!=3) 
              stop("length of starting values of parameters, 'parStart', for ggamma model must be 3.")
           good <- parStart > 0
           if (sum(!good) > 0) stop("starting values of parameters must be positive")
           else { userStart <- TRUE
                  thetaStart <- log(parStart)                                
                } 
      } else { ## if model =="lognorm"
               if (length(parStart)!=2) 
                  stop("length of starting values of parameters, 'parStart', for lognorm model must be 2.")
               if (parStart[2] < 0)
                   stop("starting value of the 2nd parameter must be positive") 
               else {  userStart <- TRUE
                       thetaStart <- parStart
                       thetaStart[2] <- log(parStart[2])  
                   } 
        } ## if model =="lognorm"     
  } ## if (!is.null(parStart))

 
  if (model=="ggamma"){
       if (sum(is.finite(lower))==3)
           lower <- log(lower)
       else if (length(lower)>1 && length(lower)<3) stop("length of 'lower' must be 3.")
       if (sum(is.finite(upper))==3)
           upper <- log(upper) 
       else if (length(upper)>1 && length(upper)<3) stop("length of 'upper' must be 3.")
  } else if (model=="lognorm"){
             if (sum(is.finite(lower))==2)
                  lower[2] <- log(lower[2])
             else if (length(lower)>1 && length(lower)<2) stop("length of 'lower' must be 2.")
             if (sum(is.finite(upper))==2)
                   upper[2] <- log(upper[2])
             else if (length(upper)>1 && length(upper)<2) stop("length of 'upper' must be 2.")
     } 

  x <- data
  if (!userStart) {## get initial values of parameters ...
      if (model=="lognorm") {
                 loglik.y.gg.mic <- loglik.y.logN.mic
                 loglik.y.grad.mic <- loglik.y.grad.logN.mic
                 thetaStart <- c(.6,log(.3))  ## rough start for lognorm
      } else thetaStart <- rep(log(2),3) ## rough start for gen gamma
      b <- optim(par=thetaStart,fn=loglik.y.gg.mic, gr=loglik.y.grad.mic, method="L-BFGS-B", x=x, lower=lower, upper=upper,...) 
      thetaStart <- b$par
  } ## if (!userStart)

  par <- c() ## initializing estimated parameters
  object <- list()
  
  if (model=="lognorm") {
         loglik.mic <- loglik.logN.mic
         dloglik.nlm.mic <- dloglik.logN.nlm.mic
  }

  if (optimizer[1]=="optim"){
           if (!(optimizer[2] %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"))){
                 warning("unknown optim() method, `Nelder-Mead' were used")
                  optimizer[2] <- "Nelder-Mead"
           }
           if (optimizer[2] %in% c("BFGS", "CG", "L-BFGS-B")){
               if (is.na(optimizer[3])){ 
                  warning("the third parameter of optimizer argument is not supplied, 
                     finite-difference approximation of the gradient were used")
                  grr <- NULL
               } else if (!(optimizer[3] %in% c("fd","grad"))){
                     warning("in the third parameter of optimizer, only `fd' and `grad'
                          options are possible, finite-difference 
                          approximation of the gradient were used")
                     grr <- NULL
               }  else if (optimizer[3] == "grad") { 
                               if (model=="lognorm") grr <- dloglik.logN.mic
                               else grr <- dloglik.mic
                       } else    grr <- NULL
           } else grr <- NULL

            b <- optim(par=thetaStart,fn=loglik.mic, gr=grr, method=optimizer[2],x=x, r=r, lower=lower, upper=upper,...)

           
           theta <- b$par
           if (model=="ggamma") {
                 par <- exp(theta)               
           } else if (model=="lognorm"){
                  par <- theta 
                  par[2] <- exp(theta[2])
             } 
           
           ll <- b$value
           iterations <- b$counts
           termcode <- b$convergence
           if (termcode == 0)
                     conv <- "Successful completion"
           else if (termcode == 1)  
                      conv <- "The iteration limit `maxit' had been reached"
           else if (termcode == 10)  
                      conv <- "Degeneracy of the Nelder-Mead simplex"
           else if (termcode == 51)  
                      conv <- "A warning from the `L-BFGS-B' method; see help for `optim' for further details"
           else if (termcode == 52)  
                      conv <- "An error from the `L-BFGS-B' method; see help for `optim' for further details"
      }   ## end if (optimizer[1]=="optim")
       else if (optimizer[1]=="nlm"){ ## nlm() with analytical derivatives...
            b <- nlm(f=dloglik.nlm.mic, p=thetaStart,iterlim=100,x=x,r=r) 
           } else if (optimizer[1]=="nlm.fd"){ ## nlm() with finite difference deriv-s...
                 b <- nlm(f=loglik.mic, p=thetaStart,iterlim=100,x=x,r=r) 
             }

   if (optimizer[1]== "nlm.fd" || optimizer[1]== "nlm") {
         theta <- b$estimate
         if (model=="ggamma")
                 par <- exp(theta)
         else {par[1] <- theta[1]; par[2] <- exp(theta[2]) }
         ll <- b$minimum
         ## dll <- b$gradient
         iterations <- b$iterations
         termcode <- b$code
         if (termcode == 1)
                conv <- "Relative gradient is close to zero, current iterate is probably solution"
         else if (termcode == 2)  
                conv <- "Successive iterates within tolerance, current iterate is probably solution"
         else if (termcode == 3)  
                conv <- "Last global step failed to locate a point lower than `estimate'. Either `estimate' is an approximate local minimum of the function or `steptol' is too small"
         else if (termcode == 4)  
                conv <- "Iteration limit exceeded"
         else if (termcode == 5)  
                conv <- "Maximum step size `stepmax' exceeded five consecutive 
                   times. Either the function is unbounded below, becomes asymptotic 
                   to a finite value from above in some direction or stepmax is too small"   
    } ## end if (optimizer[1]== "nlm.fd" || optimizer[1]== "nlm")

  object$loglik <- ll
  object$par <- par  ## estimated parameters on original scale
  if (model=="lognorm") names(object$par) <- c("mu_fibers", "sig_fibers")
  else names(object$par) <- c("b_fibers", "d_fibers", "k_fibers")
  object$logpar <- theta ## estimated transformed parameters
  
   ## getting things after fit ...
   if (model=="lognorm"){
      object$cov.logpar <- solve(-hessian.logN.mic(theta=object$logpar, x=x, r=r))
      grad <- c(1,par[2])
      object$cov.par <- diag(grad)%*%object$cov.logpar%*%diag(grad)   
   } else{
         object$cov.logpar <- solve(-hessian.mic(theta=object$logpar, x=x, r=r))
         grad <- par
         object$cov.par <- diag(grad)%*%object$cov.logpar%*%diag(grad)         
    }

  ## getting estimates of the mean lengths of fiber in a standing tree...
  if (model=="lognorm") fy.fibers <- fyj.logN
  I <- integrate(function(y,th.fibers,r) fy.fibers(y,th.fibers)/(pi*r+2*y), lower=0,upper=20, th.fibers=par,r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
  object$mu.fibers <- (1/I$value -pi*r)/2

  
  object$termcode <- termcode
  object$conv <- conv
  object$iterations <- iterations 

  object
} ## end fled.micro


##########################################################

fled <- function(data=stop("No data supplied"), data.type="ofa", r=2.5, model="ggamma", method="ML", parStart=NULL, fixed=NULL,optimizer=c("optim","L-BFGS-B","grad"), lower=-Inf, upper=Inf,cluster=1, ...){
## function to estimate wood fiber length from increment cores...
## * data - a vector of cell length from increment cores
## * data.type - type of data supplied: "ofa" (default) or "microscopy"
## * r - radius of the increment core
## * model - if model="ggamma" then the length distribution is assumed to be generalized gamma and maximum likelihood method is used for model parameter estimation; if model="lognorm" then log normal distribution is assumed and the parameters are estimated by a stochastic version of the EM algorithm
## * method - either "ML" for maximum likelihood method (default) or "SEM" for a stochastic version of the EM algorithm used only with log normal model ("SEM" is currently turned off)
## * parStart - numerical vector of starting values of parameters (or fixed values for ggamma model)
## * fixed - TRUE/FALSE vector of seven used to specify if some parameters of ggamma model should be fixed
## * optimizer - numerical optimization method used to minimize 'minus' loglik: 'optim' or 'nlm' or 'nlm.fd' (nlm based on finite-difference approximation of the derivatives). If optimizer=="optim" then the second argument specifies the numerical method to be used in 'optim' ("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN". The third element indicates whether the finite difference approximation should be used ("fd") or analytical gradient ("grad") for the "BFGS", "CG" and "L-BFGS-B" methods. 
## * lower, * upper - lower and upper bounds on parameters estimates for the "L-BFGS-B" method. these are supplied on original scale. the function then transform them to ensure positiveness and [0,1] restrictions...
## * cluster is either '0' for no parallel computing to be used; or '1' (default) for one less than
## the number of cores; or user-supplied cluster on which to do estimation; used for analyzing OFA data only.(a cluster here can be some cores of a single machine).

  if (!is.vector(data) || is.list(data)){
     data<- unlist(data); data <- as.vector(data,mode='numeric') 
  }
  if (length(data)<7) stop("not enough data supplied")
  bad <- data <= 0
  if (sum(bad)>0) {
      warning("some data are nonpositive, removed.")
      data <- data[!bad] 
  }
  bad <- data >= 2*r
  if (sum(bad)> 0){
      warning("some data values are more than the diameter, removed.")
      data <- data[!bad] 
  }
  if (is.na(model)) model <- "ggamma"
  if (!model%in%c("ggamma","lognorm")) stop("unknown fiber length distribution.")  
  if (!data.type%in%c("ofa","microscopy")) stop("unknown type of data.")
  if (data.type=="microscopy" & !is.null(fixed))
           stop("fixing parameters is possible only with 'Kajjani' data and 'ggamma' model.")
    
  object <- list()
  if (data.type=="ofa") 
     object <- fled.kajaani(data=data, r=r, model=model, method=method, parStart=parStart, fixed=fixed, optimizer=optimizer, lower=lower, upper=upper, cluster=cluster,...)
  else 
     object <- fled.micro(data=data, r=r, model=model, parStart=parStart, optimizer=optimizer, lower=lower, upper=upper,...)
 
  object$data.type <- data.type
  object$fixed <- fixed
  object$model <- model
  object$method <- method
  object$n <- length(data)
  object$data <- data ## needs for plotting histogram with plot.fled
  object$r <- r
  class(object)<-"fled"
  object
} ## end fled


#########################
## loading functions...
##########################

print.fiberLD.version <- function()
{ library(help=fiberLD)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("This is fiberLD ",version,".",sep="")
  packageStartupMessage(hello)
}


.onAttach <- function(...) { 
  print.fiberLD.version()
}



