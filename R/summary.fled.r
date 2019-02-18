###################################################################
## printing the results of the scam (clone of print.gam())...    ##
###################################################################

print.fled <- function (x,...) 
## default print function for fled objects
{
    cat("\n")
    if (x$data.type=="ofa")
       cat("Increment core data (all fiber and fine lengths in the core)\n", sep = "")
    else cat("Microscopy data (uncut fibers in the core)\n", sep = "")
    cat("\n")
    if (x$model=="ggamma")
       cat("Model: Generalized gamma\n", sep = "")
    else cat("Model: Log normal\n", sep = "")
    if (!is.null(x$fixed)) 
        cat("\nFixed model parameters:\n", x$fixed, "\n")
    cat("\nModel parameters:\n", formatC(x$par,digits = 4), "\n", sep = " ")
    cat("\n'-'Loglik = ", round(x$loglik, digits=3), "  n = ", x$n, "\n",sep="")
    cat("\n")
    invisible(x)
}



#####################
## summary ...
#####################
Dmu.fin <- function(th.fines, r,int.val, model){
## function to calculate gradient of the mean fine/fiber lengths estimates...
## * th.fines - vector of parameters on original scale
   pir <- pi*r
   D <- matrix(NA,length(th.fines), 1)
   if (model=="lognorm") {
     int <- integrate(function(y,th.fines,pir) DlogN.mu(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
     D[1] <- -int$value/2/int.val^2## derivative wrt mu_1
     int <- integrate(function(y,th.fines,pir) DlogN.logsig(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
     D[2] <- -int$value/2/int.val^2## derivative wrt log(sig1)     
  } else {
     int <- integrate(function(y,th.fines,pir) Dggamma.b(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
     D[1] <- -int$value/2/int.val^2## derivative wrt theta_2
     int <- integrate(function(y,th.fines,pir) Dggamma.d(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
     D[2] <- -int$value/2/int.val^2## derivative wrt theta_3
     int <- integrate(function(y,th.fines,pir) Dggamma.k(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
     D[3] <- -int$value/2/int.val^2## derivative wrt theta_3     
   }
  D
} ## end Dmu.fin

gradEW.k <- function(th.fines, r, int.val, model, mu.fin,k){
## function to calculate gradient of the EW^k estimate...
## * th.fines - vector of parameters on original scale
   pir <- pi*r
   d <- matrix(NA,length(th.fines), 1)
   Dmu <- Dmu.fin(th.fines, r,int.val, model)
   if (model=="lognorm") {
     d[1] <- integrate(function(y,th.fines,pir,mu.fin,dmu, k) y^k/(pir+2*y)*(2*fyj.logN(y,th.fines)*dmu+ DlogN.mu(y,th.fines)*(pir+2*mu.fin)), lower=0,upper=20, th.fines=th.fines,pir=pir, mu.fin=mu.fin, dmu=Dmu[1], k=k, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
     d[2] <- integrate(function(y,th.fines,pir,mu.fin,dmu,k) y^k/(pir+2*y)*(2*fyj.logN(y,th.fines)*dmu+ DlogN.logsig(y,th.fines)*(pir+2*mu.fin)), lower=0,upper=20, th.fines=th.fines,pir=pir, mu.fin=mu.fin, dmu=Dmu[2],k=k, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value   
  } else {
       d[1] <- integrate(function(y,th.fines,pir,mu.fin,dmu, k) y^k/(pir+2*y)*(2*fy.fines(y,th.fines)*dmu+ Dggamma.b(y,th.fines)*(pir+2*mu.fin)), lower=0,upper=20, th.fines=th.fines,pir=pir, mu.fin=mu.fin, dmu=Dmu[1], k=k,stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
       d[2] <- integrate(function(y,th.fines,pir,mu.fin,dmu,k) y^k/(pir+2*y)*(2*fy.fines(y,th.fines)*dmu+ Dggamma.d(y,th.fines)*(pir+2*mu.fin)), lower=0,upper=20, th.fines=th.fines,pir=pir, mu.fin=mu.fin, dmu=Dmu[2],k=k, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
       d[3] <- integrate(function(y,th.fines,pir,mu.fin,dmu,k) y^k/(pir+2*y)*(2*fy.fines(y,th.fines)*dmu+ Dggamma.k(y,th.fines)*(pir+2*mu.fin)), lower=0,upper=20, th.fines=th.fines,pir=pir, mu.fin=mu.fin, dmu=Dmu[3],k=k, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
    }
  d
} ## end gradEW.k
 



summary.fled <- function (object,...) 
## summary method for fled object...
{
  ## table for model parameters...
     p.table <- matrix(NA, 2,length(object$par))
     if (object$data.type=="microscopy"){
        p.table[1,] <- object$par ## para-s estimates
        p.table[2,] <- diag(object$cov.par)^.5 ## approx. std errors of par estimates
        dimnames(p.table) <- list(c("Estimate","Std. Error"),c(names(object$par)))
     } else {
         p.table[1,] <- c(object$par[-1], object$par[1]) ## para-s estimates
         std <- diag(object$cov.par)^.5 ## approx. std errors of par estimates
         if (!is.null(object$fixed)){
            std1 <-rep(0,length(object$par))
            std1[!object$fixed] <- std
            std <- std1
         }
         p.table[2,] <- c(std[-1],std[1])
         dimnames(p.table) <- list(c("Estimate","Std. Error"),c(names(object$par[-1]),names(object$par[1])))         
       }
  ## table for summary statistics for cell lengths in the increment core...
  if (object$data.type=="microscopy"){
     ss.table <- NULL
  } else {  ## for Kajanni data...
   if (object$model=="lognorm") {
     emu1 <- exp(object$par[2])
     esig1.sq <- exp(object$par[3]^2)
     emu2 <- exp(object$par[4])
     esig2.sq <- exp(object$par[5]^2)
     ss.table <- matrix(NA,7,2)
     ss.table[1,] <- c(emu1*esig1.sq^.5, emu2*esig2.sq^.5) ## fines and fibers means 
     ss.table[2,] <- c((esig1.sq-1)*emu1^2*esig1.sq, (esig2.sq-1)*emu2^2*esig2.sq)^.5 ## std. deviations of fines and fibers
     ss.table[3,] <- c((esig1.sq+2)*(esig1.sq-1)^.5, (esig2.sq+2)*(esig2.sq-1)^.5 ) ## skewness of fines and fibers
     ss.table[4,] <- c( (esig1.sq^4+2*esig1.sq^3+3*esig1.sq^2-6), (esig2.sq^4+2*esig2.sq^3+3*esig2.sq^2-6) ) ## excess kurtosis of fines and fibers
     ss.table[c(5:7),1] <- emu1*exp(object$par[3]*qnorm(c(.25,.5,.75))) ## 1st, 2nd, 3rd quartile of fines 
     ss.table[c(5:7),2] <- emu2*exp(object$par[5]*qnorm(c(.25,.5,.75))) ## 1st, 2nd, 3rd quartile of fibers 
     dimnames(ss.table) <- list(c("Mean","Std.dev.","Skewness","Ex.kurt.","1st Qu","Median","3rd Qu"), c("fines", "fibers") )
    } else { ## if model=="ggamma"
        ss.table <- matrix(NA,4,2)
        b1 <- object$par[2]; d1 <- object$par[3]; k1 <- object$par[4]
        b2 <- object$par[5]; d2 <- object$par[6]; k2 <- object$par[7]
        gak1 <- gamma(k1); gak2 <- gamma(k2)
        id1 <- 1/d1; id2 <- 1/d2
        ss.table[1,] <- c(b1*gamma(k1+id1)/gak1, b2*gamma(k2+id2)/gak2) ## fines and fibers means 
        ss.table[2,] <- c( b1^2*(gamma(k1+2*id1) - gamma(k1+id1)^2/gak1)/gak1 , 
           b2^2*(gamma(k2+2*id2) - gamma(k2+id2)^2/gak2)/gak2 )^.5 ## std. deviations of fines and fibers
        ss.table[3,1] <- (b1^3*gamma(k1+3*id1)/gak1 - 3*ss.table[1,1]*b1^2*gamma(k1+2*id1)/gak1 + 2*ss.table[1,1]^3 )/ss.table[2,1]^3 ## skewness of fines 
        ss.table[3,2] <- (b2^3*gamma(k2+3*id2)/gak2 - 3*ss.table[1,2]*b2^2*gamma(k2+2*id2)/gak2 + 2*ss.table[1,2]^3 )/ss.table[2,2]^3 ## skewness of fibers 
        ss.table[4,1] <- (b1^4*gamma(k1+4*id1)/gak1 - 4*ss.table[1,1]*b1^3*gamma(k1+3*id1)/gak1 + 6*ss.table[1,1]^2*ss.table[2,1]^2 +3*ss.table[1,1]^4  )/ss.table[2,1]^4 ## kurtosis of fines 
        ss.table[4,2] <- (b2^4*gamma(k2+4*id2)/gak2 - 4*ss.table[1,2]*b2^3*gamma(k2+3*id2)/gak2 + 6*ss.table[1,2]^2*ss.table[2,2]^2 +3*ss.table[1,2]^4  )/ss.table[2,2]^4 ## kurtosis of fibers  

        good <- (k1*d1 >1) && (k2*d2 >1)
        if (good){
          ss.table <- rbind(ss.table, c(b1*(k1-id1)^id1, b2*(k2-id2)^id2) ) ## getting modes of fines and fibers
          dimnames(ss.table) <- list(c("Mean","Std.dev.","Skewness","Kurtosis", "Mode"), c("fines", "fibers") )
        } else
           dimnames(ss.table) <- list(c("Mean","Std.dev.","Skewness","Kurtosis"), c("fines", "fibers") )
     } ## end if model=="ggamma"
  } # end else for Kajanni data  

  ## table for summary statistics for FINE lengths in a standing tree...
  par <- object$par
  pir <- pi*object$r
  if (object$model=="lognorm") {
      fy.fines <- fyj.logN
      if (object$data.type=="microscopy"){
           th.fines <- par
           ind <- c(1:2)
      } else {
           th.fines <- par[2:3]
          ind <- c(2:3)
        }
  } else if (object$model=="ggamma") { 
            if (object$data.type=="microscopy"){
               th.fines <- par
               ind <- c(1:3)
            } else { 
                th.fines <- par[2:4]
                ind <- c(2:4)
              }
    }
  fixed.fin <- rep(FALSE, length(ind))
  fixed.fb <- rep(FALSE, length(ind))
  if (!is.null(object$fixed)){
     cov.logpar <- matrix(0,length(object$par),length(object$par))
     cov.logpar[!object$fixed,!object$fixed] <- object$cov.logpar
     fixed.fin <- object$fixed[2:4]
     fixed.fb <- object$fixed[5:7]
  } else cov.logpar <- object$cov.logpar
  fixed.fin <- which(fixed.fin)
  fixed.fb <- which(fixed.fb)

  if (object$data.type=="microscopy"){
     w.fine <- NULL
  } else {  ## for Kajaani data...
    # if (object$model=="lognorm") {
    #    fy.fines <- fyj.logN
    #    th.fines <- par[2:3]
    #    ind <- c(2:3)
    # } else { th.fines <- par[2:4]
    #         ind <- c(2:4)
    #       }
     I=integrate(function(y,th.fines,pir) fy.fines(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
     mu.fin <- (1/I$value -pir)/2
     EW2.fin <- integrate(function(y,th.fines,pir,mu.fin) fy.fines(y,th.fines)*y^2*(pir+2*mu.fin)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,mu.fin=mu.fin, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
     EW3.fin <- integrate(function(y,th.fines,pir,mu.fin) fy.fines(y,th.fines)*y^3*(pir+2*mu.fin)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,mu.fin=mu.fin, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
     EW4.fin <- integrate(function(y,th.fines,pir,mu.fin) fy.fines(y,th.fines)*y^4*(pir+2*mu.fin)/(pir+2*y), lower=0,upper=20, th.fines=th.fines,pir=pir,mu.fin=mu.fin, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
     sig.fin <- (EW2.fin-mu.fin^2)^.5

     w.fine <- matrix(NA,2,4)
     w.fine[1,1:2] <- c(mu.fin,sig.fin)
     w.fine[1,3] <- (EW3.fin-3*EW2.fin*mu.fin+2*mu.fin^3)/sig.fin^3 
     w.fine[1,4] <- (EW4.fin-4*EW3.fin*mu.fin+6*mu.fin^2*sig.fin^2+3*mu.fin^4)/sig.fin^4 
     dimnames(w.fine) <- list(c("Estimate", "Std. Error"), c("Mean", "Std.dev.", "Skewness", "Kurtosis"))
     ## using delta method to calculate the std. errors of the summary statistics...
     
     Dmu <- Dmu.fin(th.fines=th.fines, r=object$r,int.val=I$value, model=object$model)
     Dmu[fixed.fin] <- 0
     w.fine[2,1] <- (t(Dmu)%*%cov.logpar[ind,ind]%*%Dmu)^.5 ## std err of the mean
     Dstd.err <- matrix(NA,length(th.fines),1)
     dEW2 <- gradEW.k(th.fines=th.fines, r=object$r,int.val=I$value, model=object$model, mu.fin=mu.fin,k=2)
     dEW2[fixed.fin] <- 0
     Dstd.err <- (dEW2-2*mu.fin*Dmu)/2/sig.fin
     w.fine[2,2] <- (t(Dstd.err)%*%cov.logpar[ind,ind]%*%Dstd.err)^.5 ## std error of sigma
     dEW3 <- gradEW.k(th.fines=th.fines, r=object$r,int.val=I$value, model=object$model, mu.fin=mu.fin,k=3)
     dEW3[fixed.fin] <- 0
     dEW4 <- gradEW.k(th.fines=th.fines, r=object$r,int.val=I$value, model=object$model, mu.fin=mu.fin,k=4)
     dEW4[fixed.fin] <- 0
     Dskew <- (dEW3-3*mu.fin*dEW2-3*EW2.fin*Dmu+6*mu.fin^2*Dmu)/sig.fin^3 - 3*Dstd.err*w.fine[1,3]/sig.fin
     w.fine[2,3] <- (t(Dskew)%*%cov.logpar[ind,ind]%*%Dskew)^.5 ## std err of skewness
     Dkurt <- (dEW4-4*EW3.fin*Dmu-4*mu.fin*dEW3+12*mu.fin*sig.fin^2*Dmu+12*mu.fin^2*sig.fin*Dstd.err +12*mu.fin^3*Dmu)/sig.fin^4 - 4*Dstd.err*w.fine[1,4]/sig.fin
     w.fine[2,4] <- (t(Dkurt)%*%cov.logpar[ind,ind]%*%Dkurt)^.5 ## std err of kurtosis
   } ## end of w.fine for Kajaani data

## checking deriv by finite differencing...
 # del <- 1e-5; theta <-th.fines; 
 # if (object$model=="ggamma")  theta <- log(th.fines)
 #  else  theta[2]<- log(th.fines[2])
 # fifmu <- rep(0,length(th.fines))  
 # findif <- rep(0,length(th.fines))
 # for (i in 1:length(th.fines)){
 #   theta1 <- theta; theta1[i] <- theta[i]+del
 #   th1 <-theta1; 
 #   if (object$model=="ggamma") th1 <- exp(theta1)               
 #    else  th1[2]<-exp(theta1[2])  
 #   I1=integrate(function(y,th.fines,pir) fy.fines(y,th.fines)/(pir+2*y), lower=0,upper=20, th.fines=th1,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.1) 
 #   mu1 <- (1/I1$value -pir)/2
 #   EW2.1<- integrate(function(y,th.fines,pir,mu.fin) fy.fines(y,th.fines)*y^2*(pir+2*mu.fin)/(pir+2*y), lower=0,upper=20, th.fines=th1,pir=pir,mu.fin=mu1, stop.on.error = F, rel.tol=.Machine$double.eps^0.1)$value
  #  sig1 <- (EW2.1-mu1^2)^.5
  #  findif[i] <- (sig1-sig.fin)/del
  #  fifmu[i] <- (mu1-mu.fin)/del
  #}
# print((Dstd.err-findif)/Dstd.err)
# print(Dstd.err); print(findif)
# print((Dmu-fifmu)/Dmu)
# print(Dmu); print(fifmu)

 ##-------------------
 ## table for summary statistics for FIBER lengths in a standing tree...
  if (object$model=="lognorm") {
      fy.fibers <- fyj.logN
      if (object$data.type=="microscopy"){
           th.fibers <- par
           ind <- c(1:2)
      } else {
           th.fibers <- par[4:5]
           ind <- c(4:5)
        }
  } else if (object$model=="ggamma") { 
            if (object$data.type=="microscopy"){
               th.fibers <- par
               ind <- c(1:3)
            } else { 
                th.fibers <- par[5:7]
                ind <- c(5:7)
              }
    }
   I=integrate(function(y,th.fibers,pir) fy.fibers(y,th.fibers)/(pir+2*y), lower=0,upper=20, th.fibers=th.fibers,pir=pir,stop.on.error = F, rel.tol=.Machine$double.eps^0.3) 
   mu.fb <- (1/I$value -pir)/2
   EW2.fb <- integrate(function(y,th.fibers,pir,mu.fb) fy.fibers(y,th.fibers)*y^2*(pir+2*mu.fb)/(pir+2*y), lower=0,upper=20, th.fibers=th.fibers,pir=pir,mu.fb=mu.fb, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   EW3.fb <- integrate(function(y,th.fibers,pir,mu.fb) fy.fibers(y,th.fibers)*y^3*(pir+2*mu.fb)/(pir+2*y), lower=0,upper=20, th.fibers=th.fibers,pir=pir,mu.fb=mu.fb, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
   EW4.fb <- integrate(function(y,th.fibers,pir,mu.fb) fy.fibers(y,th.fibers)*y^4*(pir+2*mu.fb)/(pir+2*y), lower=0,upper=20, th.fibers=th.fibers,pir=pir,mu.fb=mu.fb, stop.on.error = F, rel.tol=.Machine$double.eps^0.3)$value
  sig.fb <- (EW2.fb-mu.fb^2)^.5

    w.fiber <- matrix(NA,2,4)
    w.fiber[1,1:2] <- c(mu.fb,sig.fb)
    w.fiber[1,3] <- (EW3.fb-3*EW2.fb*mu.fb+2*mu.fb^3)/sig.fb^3 
    w.fiber[1,4] <- (EW4.fb-4*EW3.fb*mu.fb+6*mu.fb^2*sig.fb^2+3*mu.fb^4)/sig.fb^4 
    dimnames(w.fiber) <- list(c("Estimate", "Std. Error"), c("Mean", "Std.dev.", "Skewness", "Kurtosis"))
    ## using delta method to calculate the std. errors of the summary statistics..
    Dmu.fb <- Dmu.fin(th.fines=th.fibers, r=object$r,int.val=I$value, model=object$model)
    Dmu.fb[fixed.fb] <- 0
    w.fiber[2,1] <- (t(Dmu.fb)%*%cov.logpar[ind,ind]%*%Dmu.fb)^.5 ## std err of the mean
    Dstd.err <- matrix(NA,length(th.fibers),1)
    dEW2 <- gradEW.k(th.fines=th.fibers, r=object$r,int.val=I$value, model=object$model, mu.fin=mu.fb,k=2)
    dEW2[fixed.fb] <- 0
    Dstd.err <- (dEW2-2*mu.fb*Dmu.fb)/2/sig.fb
    w.fiber[2,2] <- (t(Dstd.err)%*%cov.logpar[ind,ind]%*%Dstd.err)^.5 ## std err of std.err
    dEW3 <- gradEW.k(th.fines=th.fibers, r=object$r,int.val=I$value, model=object$model, mu.fin=mu.fb,k=3)
    dEW3[fixed.fb] <- 0
    dEW4 <- gradEW.k(th.fines=th.fibers, r=object$r,int.val=I$value, model=object$model, mu.fin=mu.fb,k=4)
    dEW4[fixed.fb] <- 0
    Dskew <- (dEW3-3*mu.fb*dEW2-3*EW2.fb*Dmu.fb+6*mu.fb^2*Dmu.fb)/sig.fb^3 - 3*Dstd.err*w.fiber[1,3]/sig.fb
    w.fiber[2,3] <- (t(Dskew)%*%cov.logpar[ind,ind]%*%Dskew)^.5 ## std err of skewness
    Dkurt <- (dEW4-4*EW3.fb*Dmu.fb-4*mu.fb*dEW3+12*mu.fb*EW2.fb*Dmu.fb+6*mu.fb^2*dEW2 - 12*mu.fb^3*Dmu.fb)/sig.fb^4 - 4*Dstd.err*w.fiber[1,4]/sig.fb
    w.fiber[2,4] <- (t(Dkurt)%*%cov.logpar[ind,ind]%*%Dkurt)^.5 ## std err of kurtosis
 
  ##----------------------
  ## getting estimate of the expected value of the cell length in a standing tree, E(W)...
  Dmean.w <- function(eps, mu1, mu2, r){
   ## gets gradient of the estimate of E(W)...
   ## * mu1 - mean of fine length (mu.fin)
   ## * mu2 - mean of fiber length (mu.fb)
         pir <- pi*r 
         A <- 2*(mu2*eps +(1-eps)*mu1)+pir
         B <- 2*mu1*mu2+mu1*pir*eps +mu2*pir*(1-eps) 
         d <- rep(NA, 3)
         d[1] <- (eps^2 -eps)*(mu2-mu1)*( pir + 2*B/A )/A ## deriv wrt tau (transformed eps)
         d[2] <- (2*mu2+pir*eps)/A - 2*(1-eps)*B/A^2   ## deriv wrt mu1
         d[3] <- (2*mu1+pir*(1-eps))/A - 2*eps*B/A^2   ## deriv wrt mu2
         d
  } ## end Dmean.w

  if (object$data.type=="microscopy"){
     mean.w <- eps.tree <- se.eps.tree <- NULL
  } else { ## for Kajaani data... 
     mean.w <- object$mu.cell 
     ## getting the proportion of fines in a standing tree...
     eps.tree <- object$prop.fines  ## par[1]*(pir+2*mean.w)/(pir+2*mu.fin)
     ## getting std. error of the above proportion...
     # if (is.null(object$fixed)){
     eps <- par[1]
     dEW.all <- Dmean.w(eps=par[1], mu1=mu.fin, mu2=mu.fb, r=object$r)
     ## Dmu, Dmu.fb
     .dEW <- rep(NA,length(par)) ## derivative of E(W) wrt to all mixed model parameters 
     .dEW[1] <- dEW.all[1] ## wrt tau/theta_1
     if (object$model=="lognorm") {
         ind1 <- c(2:3); ind2 <- c(4:5)
     } else {ind1 <- c(2:4); ind2 <- c(5:7)}
     .dEW[ind1] <- dEW.all[1]*Dmu
     .dEW[ind2] <- dEW.all[2]*Dmu.fb
     deps.tree <- rep(NA,length(par)) ## gradient of the proportion of fines in a standing tree
     term <- pir+2*mu.fin
     deps.tree[1] <- (eps-eps^2)*(pir +2* mean.w)/term +
               eps*2*.dEW[1]/term
     deps.tree[ind1] <- 2*eps*(.dEW[ind1]*term -Dmu*(pir +2* mean.w))/term^2
     deps.tree[ind2] <- 2*eps*.dEW[ind2]/term
     se.eps.tree <-(t(deps.tree)%*%cov.logpar%*%deps.tree)^.5 ## std err of eps.tree
     # } else se.eps.tree <- NULL
   } ## end else for Kajaani data

  ret <- list(model=object$model, n=object$n, loglik=object$loglik, fixed=object$fixed,
    conv=object$conv, p.table=p.table, ss.table=ss.table,w.fine=w.fine, 
    w.fiber=w.fiber, mean.w=mean.w, eps.tree=eps.tree, se.eps.tree=se.eps.tree, data.type=object$data.type, method=object$method)

  class(ret)<-"summary.fled"
  ret
} ## end summary.fled


print.summary.fled <- function(x, digits = max(3, getOption("digits") - 3), ...)
## print method for fled summary method...
{  cat("\n")
   if (x$data.type=="ofa")
      cat("Increment core data (all fiber and fine lengths in the core)\n", sep = "")
   else cat("Microscopy data (uncut fibers in the core)\n", sep = "")
   cat("\n")
   if (x$model=="ggamma")
       cat("Model: Generalized gamma", "   Method: ML", "\n",sep="")
   else cat("Model: Log normal", "   Method: ", x$method, "\n",sep="")
   if (!is.null(x$fixed)) 
        cat("\nFixed model parameters:\n", x$fixed, "\n")
  # cat("\nn = ", x$n, "\n", sep = "")

   cat("\nModel parameters:\n")
   printCoefmat(x$p.table, digits = digits, na.print = "NA", ...)

   ## table for summary statistics for length of cells in the increment core...
 #  cat("\nSummary statistics for cell lengths in the increment core:\n")
 #  printCoefmat(x$ss.table, digits = digits, na.print = "NA",cs.ind=1, ...)
   ## table for summary statistics for length of fiber in a standing tree...
   cat("\nSummary statistics for FIBER lengths in the standing tree:\n")
   printCoefmat(x$w.fiber, digits = digits, na.print = "NA",cs.ind=1, ...)

   if (x$data.type=="ofa"){
     ## table for summary statistics for length of fine in a standing tree...
     cat("\nSummary statistics for FINE lengths in the standing tree:\n")
     printCoefmat(x$w.fine, digits = digits, na.print = "NA",cs.ind=1, ...)
     ## printing proportion of fines in a standing tree...
     cat("\nProportion of fines in the standing tree: ", round(x$eps.tree,digits=2),
      sep = "")
     # if (is.null(x$fixed)) 
     cat(" (Std.error = ",round(x$se.eps.tree,digits=3), ")", sep ="") 
   }
   cat("\n")
   cat("\n'-'Loglik = ", round(x$loglik,digits=3), "   Sample size: n = ", x$n, "\n",sep="")
   if (x$method=="ML")
      cat("\nConvergence: ", x$conv, "\n", sep = "")

   invisible(x)
} ## end print.summary.fled




