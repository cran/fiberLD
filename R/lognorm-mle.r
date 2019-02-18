## R routines for wood fiber length determination based on log normal distribution with MLEs ... 

# Fines/fibers Y densities...
fyj.logN <- function(x,th.j){
## gets density of fine or fiber length as density function of log normal distribution
## * x - vector of quantiles of true length of fine
## * th.j- lognormal parameters: (mu, sigma)
     exp(-(log(x)-th.j[1])^2/2/th.j[2]^2)/x/th.j[2]/(2*pi)^.5  
}

# getting density of observed FINE/FIBER length, fines/fiber X, ...
integrand.j.logN <- function(y,x,th.j,r){
## integrand needed for calculating X fine density 
  fx.given.y(x,y,r)*fyj.logN(y,th.j)
}

fxj.logN <- function(X,th.j,r){
## gets density of observed fine/fiber length (X fines/fibers)... 
  z=c()
  for(i in 1:length(X)) {
      x=X[i]
      I=integrate(integrand.j.logN, lower=x, upper=100, x=x, th.j=th.j, 
         r=r, stop.on.error = F, rel.tol=.Machine$double.eps^.3)
      z[i]=p.uc(x,r)*fyj.logN(x,th.j)+I$value
  }
  z
} 

 fxj.logN.1 <- function(x,th.j,r){
  ## gets density of a single observation of fine/fiber length...  
  ## * th.j- vector of two parameters (mu,sigma)
   I=integrate(integrand.j.logN, lower=x, upper=20, x=x, th.j=th.j, 
         r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)
   z=p.uc(x,r)*fyj.logN(x,th.j)+I$value
  } 

 
# Mixture density (with transformed parameters)...
fx.logN.1 <-function(y,theta, r,log=TRUE){
## mixture density of a single value of the observed cell, fine and fiber...
## * theta - vector of five transformed (mu's stay untransformed) parameters(eps,mu.fines,sigma.fines,mu.fibers,sigma.fibers) 
## * eps - a scalar, proportion of fines in the increment core
  theta[c(3,5)] <- exp(theta[c(3,5)])
  eps <- exp(theta[1])/(1+exp(theta[1])) 
 
  if(y>=(2*r) || y<=0) z <- 0
  else{
     z<- (1-eps)*fxj.logN.1(y,th.j=theta[4:5],r) + eps*fxj.logN.1(y, th.j=theta[2:3],r)
  }
  if (log) z <- log(z+1e-7) else z <- z+1e-7
  z
} ## fx.logN.1


 lx.logN <- function(x,theta,r, log=TRUE){
    ## vector of mixture densities of n cells length...
       sapply(x, FUN=fx.logN.1, theta=theta,r=r,log=log) ## returns a vector
 }
  
 parlx.logN <- function(cl,x,theta,r,log=TRUE){
      ## vector of mixture densities of n cells length...
       tmp <- parLapply(cl,x, fun=fx.logN.1, theta=theta,r=r,log=log) 
       unlist(tmp)
 }  

## getting log likelihood with parallel computation....
loglik.logN <- function(cl=NULL,theta,x, r){
## minus log likelihood function of log normal model (mixture density)... 
## * theta - vector of transformed parameters     
   ## (eps,mu.fines,sig.fines,mu.fibers,sig.fibers)
## * eps - a scalar, proportion of fines in the increment core
## * x - cell length
   if (is.null(cl))   
        ll <- sum(lx.logN(x=x,theta=theta, r=r,log=TRUE))
   else       
      ll <- sum(parlx.logN(cl=cl,x=x,theta=theta,r=r,log=TRUE))    
  -ll
} ## loglik.logN


DlogN.mu <- function(y,th){
    ## gets derivative of log normal density wrt mu...
    ## th - vector of two para-s
       fyj.logN(y,th)*(log(y)-th[1])/th[2]^2
    } 
DlogN.logsig <- function(y,th){
    ## gets derivative of log normal density wrt log(sig)...
       fyj.logN(y,th)*((log(y)-th[1])^2/th[2]^2 - 1)     
} 
    

Dmixure.full.logN <- function(x, theta,r){
## gets gradient of the LOG mixure density of the full model for a single observation wrt to all five transformed para-s
## * theta - vector of transformed parameters
    ## (eps,mu.fines,sig.fines,mu.fibers,sig.fibers)
    th <- theta
    theta[c(3,5)] <- exp(theta[c(3,5)])
    eps <- exp(theta[1])/(1+exp(theta[1])) 
    dmix <- rep(NA,5)
    dmix[1] <- eps^2/exp(th[1])*(fxj.logN.1(x=x,th.j=theta[2:3],r=r) -
                  fxj.logN.1(x=x,th.j=theta[4:5],r=r)) ## deriv wrt theta[1]
    
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*DlogN.mu(y,th), lower=x,upper=20, 
           x=x, th=theta[2:3],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)  
    dmix[2] <- eps*(p.uc(x,r=r)*DlogN.mu(y=x,th=theta[2:3]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*DlogN.logsig(y,th), lower=x,upper=20, 
           x=x, th=theta[2:3],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)
    dmix[3] <- eps*(p.uc(x,r=r)*DlogN.logsig(y=x,th=theta[2:3]) + I$value)
    
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*DlogN.mu(y,th), lower=x,upper=20, 
          x=x, th=theta[4:5],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)
    dmix[4] <- (1-eps)*(p.uc(x,r=r)*DlogN.mu(y=x,th=theta[4:5]) + I$value)
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*DlogN.logsig(y,th), lower=x,upper=20, 
           x=x, th=theta[4:5],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)
    dmix[5] <- (1-eps)*(p.uc(x,r=r)*DlogN.logsig(y=x,th=theta[4:5]) + I$value)
    dmix/fx.logN.1(y=x,theta=th, r=r,log=FALSE) 
} ## end Dmixure.full.logN



 grad.mix.full.logN <- function(x,theta,r){
  ## gradient matrix of mixture densities of n cells...
      sapply(x, FUN=Dmixure.full.logN, theta=theta,r=r) ## returns a matrix
 } 

 pargrad.mix.full.logN <- function(cl,x,theta,r){
  ## matrix of mixture densities of n cells length...
      tmp <- parLapply(cl,x, fun=Dmixure.full.logN, theta=theta,r=r) 
      tmp <- unlist(tmp)
      matrix(tmp,length(theta),length(x))
 }  


dloglik.logN <- function(cl,theta, x, r){
## gets the gradient of the log likelihood function of full mixure model...
## * theta - vector of transformed parameters 
   if (is.null(cl))    
       dll <- - rowSums(grad.mix.full.logN(x=x,theta=theta,r=r))  
   else      
      dll <- -rowSums(pargrad.mix.full.logN(cl=cl,x=x,theta=theta,r=r))     
   dll
} ## end dloglik.logN




dloglik.logN.nlm <- function(cl,theta, x, r){
## needed for nlm() function to get the loglik and its gradient as an attribute...
  ll <- loglik.logN(cl,theta,x,r)
  attr(ll,"gradient") <- dloglik.logN(cl,theta,x,r)
  ll
}

## checking derivatives by finite differences...
 ## checking deriv of log mixure for a single observ...
#theta <- c(0,-1.75,log(1.08), .6,log(0.28)) 
#y <- 0.01
#dmix <- Dmixure.full.logN(x=y,theta=theta,r=6)
#del <- 1e-6
#findif <- rep(0,5)
#for (i in 1:5){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (fx.logN.1(y=y,theta=theta1, r=6,log=TRUE)-fx.logN.1(y=y,theta=theta, log=TRUE))/del
#}
#(dmix-findif)/dmix
 ## checking deriv of log likelihood ...
#theta <- c(0,-1.75,log(1.08), .6,log(0.28)) 
#y <- s[1:10]
#dll <- dloglik.logN(cl=makeCluster(detectCores()),theta=theta,x=y,r=6)
#del <- 1e-5
#findif <- rep(0,7)
#for (i in 1:7){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik.logN(cl=NULL,theta=theta1,x=y, r=6)-loglik.logN(cl=NULL,theta=theta,x=y,r=6))/del
#}
#(dll-findif)/dll


 D2logN.mu <- function(y,th){
 ## gets 2nd order partial derivatives of log normal density wrt mu...
 ## th - vector of two para-s
     DlogN.mu(y,th)*(log(y)-th[1])/th[2]^2- fyj.logN(y,th)/th[2]^2    
 }    
 D2logN.logsig <- function(y,th){
 ## gets 2nd order partial derivatives of log normal density wrt log sig...
 ## th - vector of two para-s
    DlogN.logsig(y,th)*( (log(y)-th[1])^2/th[2]^2 - 3) -2*fyj.logN(y,th)    
 }    
 D2logN.musig <- function(y,th){
 ## gets cross 2nd partial deriv of log normal density wrt mu logsig...
 ## th - vector of two para-s
    DlogN.logsig(y,th)*(log(y)-th[1])/th[2]^2 - 2*DlogN.mu(y,th)     
 }    

D2mixure.full.logN <- function(x, theta,r){
## gets Hessian matrix of the LOG mixure densty of the full model for a single observation wrt to all five transformed para-s
## * theta - vector of transformed parameters
    ## (eps,mu.fines,sig.fines,mu.fibers,sig.fibers)
    th <- theta
    theta[c(3,5)] <- exp(theta[c(3,5)])
    eps <- exp(theta[1])/(1+exp(theta[1])) 

  #  D2logN <- function(y,th){
    ## gets 2nd order partial derivatives of log normal density ...
    ## th - vector of two para-s
  #     d2 <- rep(NA, 3)
  #     logth <- (log(y)-th[1])/th[2]
  #     d2[1] <- DlogN.mu(y,th)*logth/th[2]- fyj.logN(y,th)/th[2]^2 ## direct 2nd deriv wrt mu
  #     d2[2] <- DlogN.logsig(y,th)*( logth^2 - 3) -2*fyj.logN(y,th) ## direct 2nd deriv wrt logsig
  #     d2[3] <- DlogN.logsig(y,th)*logth/th[2] - 2*DlogN.mu(y,th) ## cross 2nd partial deriv   
  #     d2
  #  } 
   
    h <- matrix(NA,5,5)
    grad <- Dmixure.full.logN(x,th,r) ## getting gradient of the log mixure
    h[1,1] <- -grad[1]^2 + (1- exp(th[1]))*grad[1]/(1+exp(th[1]))
    fm <- 1/fx.logN.1(y=x,theta=th, r=r,log=FALSE)
   # epa <- eps^2/exp(th[1])
   # h[1,2] <- h[2,1] <- -grad[1]*grad[2] + epa*fm*DlogN.mu(y=x,th=theta[2:3])
   # h[1,3] <- h[3,1] <- -grad[1]*grad[3] + epa*fm*DlogN.logsig(y=x,th=theta[2:3])
   # h[1,4] <- h[4,1] <- -grad[1]*grad[4] - epa*fm*DlogN.mu(y=x,th=theta[4:5])
   # h[1,5] <- h[5,1] <- -grad[1]*grad[5] - epa*fm*DlogN.logsig(y=x,th=theta[4:5])
    ind <- c(2:3)
    h[1,ind] <- h[ind,1] <- -grad[1]*grad[ind] + grad[ind]/(1+exp(th[1]))
    ind <- c(4,5)
    h[1,ind] <- h[ind,1] <- -grad[1]*grad[ind] - eps^2/exp(th[1])*grad[ind]/(1-eps)
    
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2logN.mu(y,th), lower=x,upper=20, 
           x=x, th=theta[2:3],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)  
    h[2,2] <- -grad[2]^2 + eps*fm*(p.uc(x,r=r)*D2logN.mu(y=x,th=theta[2:3]) + I$value)

    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2logN.musig(y,th), lower=x,upper=20, 
           x=x, th=theta[2:3],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
    h[2,3] <- h[3,2] <- -grad[2]*grad[3] + eps*fm*(p.uc(x,r=r)*D2logN.musig(y=x,th=theta[2:3]) + I$value)
    ind <-c(4,5)
    h[2,ind] <- h[ind,2] <- -grad[2]*grad[ind]
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2logN.logsig(y,th), lower=x, upper=20, x=x, th=theta[2:3],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)  
    h[3,3] <- -grad[3]^2 + eps*fm*(p.uc(x,r=r)*D2logN.logsig(y=x,th=theta[2:3]) + I$value)
    ind <-c(4,5)
    h[3,ind] <- h[ind,3] <- -grad[3]*grad[ind]
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2logN.mu(y,th), lower=x, upper=20, x=x, th=theta[4:5],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)  
    h[4,4] <- -grad[4]^2 + (1-eps)*fm*(p.uc(x,r=r)*D2logN.mu(y=x,th=theta[4:5]) + I$value)
    
    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2logN.musig(y,th), lower=x,upper=20, 
           x=x, th=theta[4:5],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
    h[4,5] <- h[5,4] <- -grad[4]*grad[5] + (1-eps)*fm*(p.uc(x,r=r)*D2logN.musig(y=x,th=theta[4:5]) + I$value)

    I=integrate(function(y,x,th,r) fx.given.y(x,y,r)*D2logN.logsig(y,th), lower=x, upper=20, x=x, th=theta[4:5],r=r,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)  
    h[5,5] <- -grad[5]^2 + (1-eps)*fm*(p.uc(x,r=r)*D2logN.logsig(y=x,th=theta[4:5]) + I$value)
    h
} ## end D2mixure.full.logN


## checking Hessian by finite differences...
 ## checking the Hessian of log mixure for a single observ...
#theta <- c(0,-1.75,log(1.08), .6,log(0.28)) 
#y <- 0.01
#hh <- D2mixure.full.logN(x=y,theta=theta,r=6)
#del <- 1e-5
#findif <- matrix(NA,5,5)
#for (i in 1:5) {
#   theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i,] <- (Dmixure.full.logN(x=y,theta=theta1, r=6) - Dmixure.full.logN(x=y,theta=theta,r=6 ))/del
#}
#(hh-findif)/hh

 hess.mix.full.logN <- function(x,theta,r){
  ## Hessian matrix of mixture densities of n cells...
      sapply(x, FUN=D2mixure.full.logN, theta=theta,r=r) ## returns a matrix
 } 

 parhess.mix.full.logN <- function(cl,x,theta,r){
  ## matrix of mixture densities of n cells length...
      tmp <- parSapply(cl,x, FUN=D2mixure.full.logN, theta=theta,r=r) 
      tmp 
 }  

hessian.logN <- function(cl,theta, x, r){
## gets the Hessian of the log likelihood function of full mixure model...
## * theta - vector of transformed parameters 
   if (is.null(cl))    
       dll <- matrix(rowSums(hess.mix.full.logN(x=x,theta=theta,r=r)), 5,5)  
   else     
      dll <- matrix(rowSums(parhess.mix.full.logN(cl=cl,x=x,theta=theta,r=r)),5,5)         
   dll
} ## end hessian.logN







#################################################################
## getting log likelihood function for `Easy life',
## assuming uncut cells (log normal distribution on X's)
## needed for parameters initialization....

fy.gg.logN<-function(x,theta, log=TRUE){
## mixture density of cell length with lognormal on X's...
## * theta - vector of transformed parameters     
  theta[c(3,5)] <- exp(theta[c(3,5)])
  eps <- exp(theta[1])/(1+exp(theta[1])) 
  z=(1-eps)*fyj.logN(x = x,th.j = theta[4:5]) +
    eps*fyj.logN(x = x,th.j = theta[2:3])
  if (log) z <- log(z+1e-7) else z <- z+1e-7
  z
}
  
  fy.logN <- function(x,theta, log=TRUE){
  ## mixture densities of n cells length...
     sapply(x, FUN=fy.gg.logN, theta=theta,log=log) ## returns a vector
  }  

  parfy.logN <- function(cl,x,theta, log=TRUE){
  ## mixture densities of n cells length...
     tmp <- parLapply(cl,x, fun=fy.gg.logN, theta=theta,log=log) ## returns a vector
     unlist(tmp)
  }  

loglik.y.logN <- function(cl,theta,x){
## minus log likelihood function of log normal model (mixture density) on X's... 
## * theta - vector of transformed parameters     
   ## (eps,mu.fines,sig.fines,mu.fibers,sig.fibers)
## * eps - a scalar, proportion of fines in the increment core
## * x - cell length  
   if (is.null(cl))   
        ll <- sum(fy.logN(x=x,theta=theta, log=TRUE))
   else       
      ll <- sum(parfy.logN(cl=cl,x=x,theta=theta,log=TRUE))    
  -ll   
} ## loglik.y.gg


 DlogN <- function(y,th){
  ## gets gradient of the gen gamma density function wrt mu and log(sigma) parameters...
  ## * th - a vector of two untransformed parameters of log normal distribution
    dln <- c(NA,NA)
    dln[1] <- fyj.logN(y,th)*(log(y)-th[1])/th[2]^2 ## wrt mu
    dln[2] <- fyj.logN(y,th)*((log(y)-th[1])^2/th[2]^2 - 1) ## wrt log(sigma)
    dln 
  }

  Dmixure1.logN <- function(y,theta){
  ## gets gradient of the LOG mixure distribution for a single observation wrt to all five transformed para-s
  ## * theta - vector of transformed parameters
    th <- theta
    theta[c(3,5)] <- exp(theta[c(3,5)])
    eps <- exp(theta[1])/(1+exp(theta[1])) 
    dmix <- rep(NA,5)
    dmix[1] <- eps^2/exp(theta[1])*(fyj.logN(x=y,th.j=theta[2:3]) - fyj.logN(x=y,th.j = theta[4:5])) ## deriv w.r.t. theta1 (logit eps)
    dmix[2:3] <- eps*DlogN(y=y, th=theta[2:3]) ## deriv w.r.t. mu.fines, log sig.fines
    dmix[4:5] <- (1-eps)*DlogN(y=y, th=theta[4:5]) ## deriv w.r.t. mu.fibers, log sig.fibers
   # dmix
    dmix/fy.gg.logN(x=y,theta=th, log=FALSE)
  }

 grad.mix.logN <- function(x,theta){
  ## matrix of gradient of mixture densities of n cells...
     sapply(x, FUN=Dmixure1.logN, theta=theta) ## returns a matrix
  }  

  pargrad.mix.logN <- function(cl,x,theta){
  ## matrix of gradient of mixture densities of n cells...
     tmp <- parLapply(cl,x, fun=Dmixure1.logN, theta=theta) 
     tmp <- unlist(tmp)
     matrix(tmp,length(theta),length(x))  ## returns a matrix
  }  


loglik.y.grad.logN <- function(cl,theta,x){
## calculates gradient of the log lik function of log normal mixure density on X's...
## * theta - vector of transformed parameters  
  if (is.null(cl))    
       dll <- -rowSums(grad.mix.logN(x=x,theta=theta))  
  else      
      dll <- -rowSums(pargrad.mix.logN(cl=cl,x=x,theta=theta))
  dll  
}


## The end of `Easy life'
#################################

####################################################
## checking derivatives by finite differences ...
 
## checking deriv of log likelihood ...
#theta <- c(0,-1.75,log(1.08), .6,log(0.28)) 
#y <- s[1:10]
#dll <- loglik.y.grad.logN(cl=NULL,theta=theta,x=y)
#del <- 1e-5
#findif <- rep(0,5)
#for (i in 1:5){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik.y.logN(cl=NULL,theta=theta1,x=y)-loglik.y.logN(cl=NULL,theta=theta,x=y))/del
#}
#(dll-findif)/dll


## end checking derivatives
##############################################

