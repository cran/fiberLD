## some density functions to be exported ....

## for Kajaani data...
dx.fibers <- function(x, par, r, model="ggamma"){
## gets density function of the fiber length observed in the increment core...
## * par - vector of fiber parameters
  if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      d <- fx.fibers(X=x,th.fibers=par,r=r)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       d <- fxj.logN(X=x,th.j=par,r=r)
   }
  d
}

dy.fibers <- function(x, par, model="ggamma"){
## gets density function of the true fiber length that at least partially appears in the increment core...
   if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      d <- fy.fibers(x=x,th.fibers=par)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       d <- fyj.logN(x=x,th.j=par)
   }
  d
}
 
dw.fibers <- function(x, par, r, model="ggamma"){  ## this function is for users convenience... 
## gets density function of the fiber length in a standing tree...
   ## getting estimates of the mean lengths of fiber in a standing tree...
   pir <- pi*r
   if (model=="lognorm") fy.fibers <- fyj.logN      
   I <- integrate(function(y,th.fibers,pir) fy.fibers(y,th.fibers)/(pir+2*y), lower=0,upper=20,  th.fibers=par,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
   mu.fibers <- (1/I$value -pir)/2
   if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fy.fibers(x=x, th.fibers=par)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fyj.logN(x=x,th.j=par)
   }
  d
}


dw.fibers.withMu <- function(x, par, r, mu.fibers, model="ggamma"){ ## this finction is used in plot.fled()...
## gets density function of the fiber length in a standing tree...
## * mu.fibers - expected value of the fiber length in a standing tree (an estimate is taken from fled() object)
   if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fy.fibers(x=x, th.fibers=par)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fyj.logN(x=x,th.j=par)
   }
  d
}




dx.mixture <- function(x, par, r, model="ggamma"){
## gets the mixture density of the cell length in the increment core...
## * par - vector of model parameters
    if (model=="ggamma"){
      if (length(par)!=7) stop("the length of 'par' vector must be seven for generalized gamma model")
      theta <- par
      theta[1] <- log(par[1]/(1-par[1]))
      theta[2:7] <- log(par[2:7])
      d <- lx(x=x,theta=theta,r=r, log=FALSE)
  }else {
       if (length(par)!=5) stop("the length of 'par' vector must be five for log normal model")
       theta <- par
       theta[1] <- log(par[1]/(1-par[1]))
       theta[c(3,5)] <- log(par[c(3,5)])  
       d <- lx.logN(x=x,theta=theta,r=r, log=FALSE)
   }
  unname(d)
}



dy.mixture <- function(x, par, model="ggamma"){
## gets the mixture density of the true length of a cell at least partially appeared in the increment core...
## * par - vector of model parameters
  if (model=="ggamma"){
      if (length(par)!=7) stop("the length of 'par' vector must be seven for generalized gamma model")
      theta <- par
      theta[1] <- log(par[1]/(1-par[1]))
      theta[2:7] <- log(par[2:7])
      d <- fy(x=x,theta=theta, log=FALSE)
  }else {
       if (length(par)!=5) stop("the length of 'par' vector must be five for log normal model")
       theta <- par
       theta[1] <- log(par[1]/(1-par[1]))
       theta[c(3,5)] <- log(par[c(3,5)])  
       d <- fy.logN(x=x,theta=theta, log=FALSE)
   }
  d
}


dw.mixture <- function(x, par, r, model="ggamma"){
## this function is for users' convenience ...
## gets the mixture density of the length of a cell in a standing tree...
## * par - vector of all model parameters
  ## getting estimates of the mean lengths of fiber and fine in a standing tree...
  if (model=="ggamma"){
      if (length(par)!=7) stop("the length of 'par' vector must be seven for generalized gamma model")       
  }else {
       if (length(par)!=5) stop("the length of 'par' vector must be five for log normal model")       
   }
  if (model=="lognorm") {
      fy.fines <- fyj.logN
      th.fines <- par[2:3]
  } else  th.fines <- par[2:4]
  pir <- pi*r
  I=integrate(function(y,th.fines,pir) fy.fines(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
  mu.fines <- (1/I$value -pir)/2
  if (model=="lognorm") {
      fy.fibers <- fyj.logN
      th.fibers <- par[4:5]
  } else  th.fibers <- par[5:7]
  I=integrate(function(y,th.fibers,pir) fy.fibers(y,th.fibers)/(pir+2*y), lower=0,upper=20, th.fibers=th.fibers,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
  mu.fibers <- (1/I$value -pir)/2
  ## getting estimate of the expected value of the cell length in a standing tree, E(W)...
  mu.cell <- (2*mu.fines*mu.fibers+par[1]*pir*mu.fines+(1-par[1])*pir*mu.fibers)/(2*(par[1]*mu.fibers+(1-par[1])*mu.fines)+pir)   
  ## getting the proportion of fines in a standing tree...
  prop.fines <- par[1]*(pir+2*mu.cell)/(pir+2*mu.fines)

    if (model=="ggamma")
        d <- prop.fines*dw.fibers.withMu(x=x, par=par[2:4], r=r, mu.fibers=mu.fines, model="ggamma") + (1-prop.fines)*dw.fibers.withMu(x=x, par=par[5:7], r=r, mu.fibers=mu.fibers, model="ggamma")
    else 
        d <- prop.fines*dw.fibers.withMu(x=x, par=par[2:3], r=r, mu.fibers=mu.fines, model="lognorm") + (1-prop.fines)*dw.fibers.withMu(x=x, par=par[4:5], r=r, mu.fibers=mu.fibers, model="lognorm")
  
  d
}


dw.mixture.withMu <- function(x, par, r, mu.fibers, mu.fines, prop.fines, model="ggamma"){
##  this finction is used in plot.fled()...
## gets the mixture density of the length of a cell in a standing tree...
## * par - vector of model parameters
## * mu.fibers, * mu.fines - expected values of the fiber/fine length in a standing tree (an estimate is taken from fled() object)
## * prop.fines - proportion of fines in a standing tree (an estimate is taken from fled() object)
    if (model=="ggamma"){
      if (length(par)!=7) stop("the length of 'par' vector must be seven for generalized gamma model")
       d <- prop.fines*dw.fibers.withMu(x=x, par=par[2:4], r=r, mu.fibers=mu.fines, model="ggamma") + (1-prop.fines)*dw.fibers.withMu(x=x, par=par[5:7], r=r, mu.fibers=mu.fibers, model="ggamma")
  }else {
       if (length(par)!=5) stop("the length of 'par' vector must be five for log normal model")
       d <- prop.fines*dw.fibers.withMu(x=x, par=par[2:3], r=r, mu.fibers=mu.fines, model="lognorm") + (1-prop.fines)*dw.fibers.withMu(x=x, par=par[4:5], r=r, mu.fibers=mu.fibers, model="lognorm")
   }
  d
}




## some densities for microscopy data ...

dw.fibers.micro.withMu <- function(x, par, r, mu.fibers, model="ggamma"){ ##  this finction is used in plot.fled()...
## gets density function of the fiber length in a standing tree for microscopy data...
## * mu.fibers - expected value of the fiber length in a standing tree (an estimate is taken from fled() object)
   if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fy.fibers(x=x, th.fibers=par)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fyj.logN(x=x,th.j=par)
   }
  d
}


dw.fibers.micro <- function(x, par, r, model="ggamma"){ ## this function is for users' convenience ...
## gets density function of the fiber length in a standing tree for microscopy data...
   ## getting estimates of the mean lengths of fiber in a standing tree...
   pir <- pi*r
   if (model=="lognorm") fy.fibers <- fyj.logN      
   I <- integrate(function(y,th.fibers,pir) fy.fibers(y,th.fibers)/(pir+2*y), lower=0,upper=20,  th.fibers=par,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
   mu.fibers <- (1/I$value -pir)/2
   if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fy.fibers(x=x, th.fibers=par)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       d <- (pi*r +2*mu.fibers)/(pi*r+2*x)*fyj.logN(x=x,th.j=par)
   }
  d
}


dx.fibers.micro <- function(x, par, r, model="ggamma"){
## gets density function of the fiber length for observed microscopy data in the increment core...
## * par - vector of fiber parameters
  if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      par <- log(par)
      d <- fmic.1(x=x,theta=par,r=r, log=FALSE)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       par[2] <- log(par[2])
       d <- f.logN.mic.1(x=x,theta=par,r=r, log=FALSE)
   }
  d
}


dy.fibers.micro <- function(x, par, model="ggamma"){
## gets density function of the true fiber length for microscopy data that at least partially appears in the increment core...
   if (model=="ggamma"){
      if (length(par)!=3) stop("the length of 'par' vector must be three for generalized gamma model")
      d <- fy.fibers(x=x,th.fibers=par)
  }else {
       if (length(par)!=2) stop("the length of 'par' vector must be two for log normal model")
       d <- fyj.logN(x=x,th.j=par)
   }
  d
}





