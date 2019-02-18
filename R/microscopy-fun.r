## R routines for wood fiber length determination using MICROSCOPY data ... 

## log normal and generalized gamma distributions with MLEs... 


#################################################################
## getting log likelihood function for `Easy life', assuming uncut cells 
## needed for parameters initialization....

## log normal distribution on X's...
lpdf.logN <- function(x,th){
## gets log density of fiber length as log density function of log normal distribution
## * x - vector of quantiles of true length of fine
## * th- lognormal parameters: (mu, sigma)
     -(log(x)-th[1])^2/2/th[2]^2 - log(x) -log(th[2]) -log(2*pi)/2  
}


loglik.y.logN.mic <- function(theta,x){
## minus log likelihood function of log normal fiber distribution on X's... 
## * theta - vector of transformed parameters, (mu.fibers,log(sig.fibers))
  theta[2] <- exp(theta[2])
  ff <- lpdf.logN(x=x,th=theta)
  -sum(ff)
} 

DloglogN <- function(x,theta){
  ## gets gradient of the LOG lognormal density function wrt mu and log(sigma) parameters...
  ## * theta - a vector of two transformed parameters of log normal distribution
    th <- theta
    th[2] <- exp(theta[2])
    dln <- c(NA,NA)
    dln[1] <- (log(x)-th[1])/th[2]^2 ## wrt mu
    dln[2] <- ((log(x)-th[1])^2/th[2]^2 - 1) ## wrt log(sigma)
    dln
}

 
loglik.y.grad.logN.mic <- function(theta,x){
## calculates gradient of the log lik function of log normal density on X's...
## * theta - vector of transformed parameters  
  dll <- sapply(x, FUN=DloglogN, theta=theta) ## returns a matrix
  -rowSums(dll)      
}


## generalized gamma distribution on X's...
loglik.y.gg.mic <- function(theta,x){
## minus log likelihood function of generalized gamma fiber distribution on X's... 
## * theta - vector of transformed parameters, (log(b),log(d),log(k))
## * eps - a scalar, proportion of fines in the increment core
  theta <- exp(theta)
  ff <- dgengamma.stacy(x=x,scale=theta[1],d=theta[2],k=theta[3], log=TRUE)
  -sum(ff)
} 

Dlog.gg <- function(x,theta){
  ## gets gradient of the LOG generalized gamma density function wrt log parameters...
  ## * theta - a vector of two transformed parameters of log normal distribution
    b <- exp(theta[1]); d <- exp(theta[2]); k <- exp(theta[3])  ## untransform
    dln <- rep(NA,3)
    xbd <- (x/b)^d
    dln[1] <- d*(xbd-k)  ## wrt log(b)
    dln[2] <- 1+log(xbd)*(k-xbd)  ## wrt log(d)
    dln[3] <- k*(log(xbd)-digamma(k))  ## wrt log(k)
    dln
}

loglik.y.grad.mic <- function(theta,x){
## calculates gradient of the log lik function of gen gamma density on X's...
## * theta - vector of transformed parameters  
  dll <- sapply(x, FUN=Dlog.gg, theta=theta) ## returns a matrix
  -rowSums(dll)      
}


## The end of `Easy life'
#################################
## checking deriv of log density for a single observ...
  ## generalized gamma density....
#theta <- c(-1,-2,-.5)
#y <- 2.0
#dmix <- Dlog.gg(x=y,theta=theta)
#del <- 1e-6
#findif <- rep(0,3)
#for (i in 1:3){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (-loglik.y.gg.mic(theta=theta1,x=y)+loglik.y.gg.mic(theta=theta,x=y))/del
#}
#(dmix-findif)/dmix
## log normal density....
#theta <- c(-1,-.5)
#y <- 2.0
#dmix <- DloglogN(x=y,theta=theta)
#del <- 1e-6
#findif <- rep(0,2)
#for (i in 1:2){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (-loglik.y.logN.mic(theta=theta1,x=y)+loglik.y.logN.mic(theta=theta,x=y))/del
#}
#(dmix-findif)/dmix
  ## checking deriv of log likelihood ...
  ## generalized gamma loglik....
#theta <- rep(log(2),3) 
#y <- c(2,3.2,1.8,4.2,3.9,4.9,3.7,2.7,5.4,2.9)
#dll <- loglik.y.grad.mic(theta=theta,x=y)
#del <- 1e-5
#findif <- rep(0,3)
#for (i in 1:3){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik.y.gg.mic(theta=theta1,x=y)-loglik.y.gg.mic(theta=theta,x=y))/del
#}
#(dll-findif)/dll
 ## log normal loglik....
#theta <- c(-2,log(2)) 
#y <- c(2,3.2,1.8,4.2,3.9,4.9,3.7,2.7,5.4,2.9)
#dll <- loglik.y.grad.logN.mic(theta=theta,x=y)
#del <- 1e-5
#findif <- rep(0,2)
#for (i in 1:2){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik.y.logN.mic(theta=theta1,x=y)-loglik.y.logN.mic(theta=theta,x=y))/del
#}
#(dll-findif)/dll



##########################################

fmic.1 <- function(x,theta,r, log=TRUE){
## log density function of a single microscopy observation based on generallized gamma model... 
## * theta - vector of transformed parameters log(b.fibers,d.fibers,k.fibers)
## * x - fiber length
 th.fibers <- exp(theta)
 Int <- integrate(function(y,x,th.fibers,r) p.uc(y,r)*fy.fibers(y,th.fibers), lower=0, upper=2*r, x=x, th.fibers=th.fibers, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)
 z <- fy.fibers(x=x,th.fibers=th.fibers)*p.uc(y=x,r=r)/Int$value
 if (log) z <- log(z+1e-7) else z <- z+1e-7
 z
}

loglik.mic <- function(theta,x,r, log=TRUE){
## minus log likelihood function of observed microscopy sample based on generallized gamma model... 
## * theta - vector of transformed parameters log(b.fibers,d.fibers,k.fibers)
## * x - fiber length
  temp <- sapply(x, FUN=fmic.1, theta=theta,r=r,log=log) ## returns a vector
  -sum(temp)
}

f.logN.mic.1 <- function(x,theta,r, log=TRUE){
## log density function of a single microscopy observation based on LOG NORMAL model... 
## * theta - vector of transformed parameters (mu.fibers,log(sig.fibers))
## * x - fiber length
 th.fibers <- theta
 th.fibers[2] <- exp(theta[2])
 Int <- integrate(function(y,x,th.fibers,r) p.uc(y,r)*fyj.logN(y,th.fibers), lower=0, upper=2*r, x=x, th.fibers=th.fibers, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)
 z <- fyj.logN(x=x,th.j=th.fibers)*p.uc(y=x,r=r)/Int$value
 if (log) z <- log(z+1e-7) else z <- z+1e-7
 z
}


loglik.logN.mic <- function(theta,x,r, log=TRUE){
## minus log likelihood function of observed microscopy sample based on LOG NORMAL model... 
## * theta - vector of transformed parameters log(b.fibers,d.fibers,k.fibers)
## * x - fiber length
  temp <- sapply(x, FUN=f.logN.mic.1, theta=theta,r=r,log=log) ## returns a vector
  -sum(temp)
}



Dlogf.mic <- function(x, theta,r){
## gets gradient of the LOG density of a single microscopy observation based on generallized gamma model wrt three log parameters... 
## * theta - vector of transformed parameters log(b.fibers,d.fibers,k.fibers)
   th <- exp(theta)
   dd <- rep(NA,3)
   Int <- integrate(function(y,x,th,r) p.uc(y,r)*fy.fibers(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   dInt <- integrate(function(y,x,th,r) p.uc(y,r)*Dggamma.b(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   dd[1] <- p.uc(y=x,r=r)*(Dggamma.b(y=x,th=th)/Int - fy.fibers(x=x,th.fibers=th)*dInt/Int^2)

   dInt <- integrate(function(y,x,th,r) p.uc(y,r)*Dggamma.d(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   dd[2] <- p.uc(y=x,r=r)*(Dggamma.d(y=x,th=th)/Int - fy.fibers(x=x,th.fibers=th)*dInt/Int^2)

   dInt <- integrate(function(y,x,th,r) p.uc(y,r)*Dggamma.k(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   dd[3] <- p.uc(y=x,r=r)*(Dggamma.k(y=x,th=th)/Int - fy.fibers(x=x,th.fibers=th)*dInt/Int^2)
   dd/fmic.1(x=x,theta=theta, r=r,log=FALSE)  
}



dloglik.mic <- function(theta, x, r){
## gets the gradient of the log likelihood function ...
## * theta - vector of transformed parameters
  gr <- sapply(x, FUN=Dlogf.mic, theta=theta,r=r) ## returns a matrix
  -rowSums(gr) 
}


Dlogf.logN.mic <- function(x, theta,r){
## gets gradient of the LOG density of a single microscopy observation based on log normal model wrt two parameters... 
## * theta - vector of transformed parameters (mu.fibers,log(sig.fibers))
   th <- theta
   th[2] <- exp(theta[2])
   dd <- rep(NA,2)
   Int <- integrate(function(y,x,th,r) p.uc(y,r)*fyj.logN(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   dInt <- integrate(function(y,x,th,r) p.uc(y,r)*DlogN.mu(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   dd[1] <- p.uc(y=x,r=r)*(DlogN.mu(y=x,th=th)/Int - fyj.logN(x,th)*dInt/Int^2)

   dInt <- integrate(function(y,x,th,r) p.uc(y,r)*DlogN.logsig(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   dd[2] <- p.uc(y=x,r=r)*(DlogN.logsig(y=x,th=th)/Int - fyj.logN(x,th)*dInt/Int^2)
   dd/f.logN.mic.1(x=x,theta=theta, r=r,log=FALSE)  
}


dloglik.logN.mic <- function(theta, x, r){
## gets the gradient of the log likelihood function for log normal model...
## * theta - vector of transformed parameters
  gr <- sapply(x, FUN=Dlogf.logN.mic, theta=theta,r=r) ## returns a matrix
  -rowSums(gr) 
}


dloglik.nlm.mic <- function(theta, x, r){
## needed for nlm() function to get the loglik and its gradient as an attribute...
  ll <- loglik.mic(theta,x,r)
  attr(ll,"gradient") <- dloglik.mic(theta,x,r)
  ll
}


dloglik.logN.nlm.mic <- function(theta, x, r){
## needed for nlm() function to get the loglik and its gradient as an attribute...
  ll <- loglik.logN.mic(theta,x,r)
  attr(ll,"gradient") <- dloglik.logN.mic(theta,x,r)
  ll
}

## checking derivatives by finite differences...
 ## checking deriv of log density for a single observ...
 ## generalized gamma distribution...
#theta <- c(-1,-2,-.5)
#y <- 2.0
#dmix <- Dlogf.mic(x=y,theta=theta,r=6)
#del <- 1e-5
#findif <- rep(0,3)
#for (i in 1:3){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (fmic.1(x=y,theta=theta1, r=6,log=TRUE)-fmic.1(x=y,theta=theta,r=6, log=TRUE))/del
#}
#(dmix-findif)/dmix

## log normal distribution...
#theta <- c(-1,-.5)
#y <- 2.0
#dmix <- Dlogf.logN.mic(x=y,theta=theta,r=6)
#del <- 1e-5
#findif <- rep(0,2)
#for (i in 1:2){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (f.logN.mic.1(x=y,theta=theta1, r=6,log=TRUE)-f.logN.mic.1(x=y,theta=theta,r=6, log=TRUE))/del
#}
#(dmix-findif)/dmix

 ## checking deriv of log likelihood ...
## generalized gamma distribution...
#theta <- rep(log(2),3) 
#y <- c(2,3.2,1.8,4.2,3.9,4.9,3.7,2.7,5.4,2.9)
#dll <- dloglik.mic(theta=theta,x=y,r=6)
#del <- 1e-5
#findif <- rep(0,3)
#for (i in 1:3){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik.mic(theta=theta1,x=y, r=6)-loglik.mic(theta=theta,x=y,r=6))/del
#}
#(dll-findif)/dll
## log normal distribution...
#theta <- c(2,log(2)) 
#y <- c(2,3.2,1.8,4.2,3.9,4.9,3.7,2.7,5.4,2.9)
#dll <- dloglik.logN.mic(theta=theta,x=y,r=6)
#del <- 1e-5
#findif <- rep(0,2)
#for (i in 1:2){
#  theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i] <- (loglik.logN.mic(theta=theta1,x=y, r=6)-loglik.logN.mic(theta=theta,x=y,r=6))/del
#}
#(dll-findif)/dll

## end checking deriv of log likelihood
##########################################

D2ll.mic.logN <- function(x, theta,r){
## gets Hessian matrix of the LOG density of the observed microscopy data for a single observation wrt to two transformed para-s for log normal model...
## * theta - vector of transformed parameters (mu.fibers,log(sig.fibers))
    th <- theta
    th[2] <- exp(th[2])
    h <- matrix(NA,2,2)
    dmu <- DlogN.mu(y=x,th=th)
    dlogsig <- DlogN.logsig(y=x, th=th)
    fm <- 1/fyj.logN(x=x,th.j=th)
    Int <- integrate(function(y,x,th,r) p.uc(y,r)*fyj.logN(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    dInt.mu <- integrate(function(y,x,th,r) p.uc(y,r)*DlogN.mu(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    dInt.logsig <- integrate(function(y,x,th,r) p.uc(y,r)*DlogN.logsig(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2logN.mu(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
     h[1,1] <- -fm^2*dmu^2 +fm*D2logN.mu(y=x, th=th) + dInt.mu^2/Int^2 - d2Int/Int
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2logN.logsig(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    h[2,2] <- -fm^2*dlogsig^2 +fm*D2logN.logsig(y=x, th=th) + dInt.logsig^2/Int^2 - d2Int/Int
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2logN.musig(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    h[1,2] <- h[2,1] <- -fm^2*dmu*dlogsig +fm*D2logN.musig(y=x, th=th) + dInt.mu*dInt.logsig/Int^2 - d2Int/Int
    h
} ## end D2ll.mic.logN


D2ll.mic <- function(x, theta,r){
## gets Hessian matrix of the LOG density of the observed microscopy data for a single observation wrt to two transformed para-s for generalized gamma model ...
## * theta - vector of transformed parameters log(b.fibers,d.fibers, k.fibers)
    th <- exp(theta)
    h <- matrix(NA,3,3)
    db <- Dggamma.b(y=x,th=th)
    dd <- Dggamma.d(y=x,th=th)
    dk <- Dggamma.k(y=x,th=th)
    fm <- 1/fy.fibers(x=x,th.fibers=th)
    Int <- integrate(function(y,x,th,r) p.uc(y,r)*fy.fibers(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    dInt.b <- integrate(function(y,x,th,r) p.uc(y,r)*Dggamma.b(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    dInt.d <- integrate(function(y,x,th,r) p.uc(y,r)*Dggamma.d(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    dInt.k <- integrate(function(y,x,th,r) p.uc(y,r)*Dggamma.k(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2ggamma.b2(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    h[1,1] <- -fm^2*db^2 +fm*D2ggamma.b2(y=x, th=th) + dInt.b^2/Int^2 - d2Int/Int
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2ggamma.bd(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value 
    h[1,2] <- h[2,1] <- -fm^2*db*dd +fm*D2ggamma.bd(y=x, th=th) + dInt.b*dInt.d/Int^2 - d2Int/Int
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2ggamma.d2(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    h[2,2] <- -fm^2*dd^2 +fm*D2ggamma.d2(y=x, th=th) + dInt.d^2/Int^2 - d2Int/Int
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2ggamma.bk(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value 
    h[1,3] <- h[3,1] <- -fm^2*db*dk +fm*D2ggamma.bk(y=x, th=th) + dInt.b*dInt.k/Int^2 - d2Int/Int
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2ggamma.dk(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value 
    h[2,3] <- h[3,2] <- -fm^2*dd*dk +fm*D2ggamma.dk(y=x, th=th) + dInt.k*dInt.d/Int^2 - d2Int/Int
    d2Int <- integrate(function(y,x,th,r) p.uc(y,r)*D2ggamma.k2(y,th), lower=0, upper=2*r, x=x, th=th, r=r, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    h[3,3] <- -fm^2*dk^2 +fm*D2ggamma.k2(y=x, th=th) + dInt.k^2/Int^2 - d2Int/Int
    h
} ## end D2ll.mic


## checking Hessian by finite differences...
## checking the Hessian of log density for a single microscopy observ for generalized gamma model...
#theta <- c(-1,-2,-.5)
#y <- 2.0
#hh <- D2ll.mic(x=y,theta=theta,r=6)
#del <- 1e-3
#findif <- matrix(NA,3,3)
#for (i in 1:3) {
#   theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i,] <- (c(Dlogf.mic(x=y,theta=theta1, r=6)) - c(Dlogf.mic(x=y,theta=theta,r=6 )))/del
#}
#(hh-findif)/hh

## checking the Hessian of log density for a single microscopy observ for log normal model...
#theta <- c(-1,-.5)
#y <- 2.0
#hh <- D2ll.mic.logN(x=y,theta=theta,r=6)
#del <- 1e-5
#findif <- matrix(NA,2,2)
#for (i in 1:2) {
#   theta1 <- theta; theta1[i] <- theta[i]+del
#  findif[i,] <- (c(Dlogf.logN.mic(x=y,theta=theta1, r=6)) - c(Dlogf.logN.mic(x=y,theta=theta,r=6 )))/del
#}
#(hh-findif)/hh

## end checking hessian
#########################

hessian.logN.mic <- function(theta, x, r){
## gets the Hessian of the log likelihood function of observed microscopy data based on log normal distribution...
## * theta - vector of transformed parameters 
  h <- sapply(x, FUN=D2ll.mic.logN, theta=theta,r=r) ## returns a matrix
  matrix(rowSums(h), 2,2)  
}


hessian.mic <- function(theta, x, r){
## gets the Hessian of the log likelihood function of observed microscopy data based on generalized gamma distribution...
## * theta - vector of transformed parameters 
  h <- sapply(x, FUN=D2ll.mic, theta=theta,r=r) ## returns a matrix
  matrix(rowSums(h), 3,3)  
}












