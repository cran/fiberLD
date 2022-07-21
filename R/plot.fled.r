## plot function ...

plot.fled <- function(x,select=NULL, density.scale="tree", rvec=NULL, xlab=NULL,ylab=NULL,  main=NULL, col=4, lwd=2, ...)   ##  ylim=NULL,xlim=NULL,          
## creates plots for the mixure model and two separate plots of fiber and fine densities...
## * x is a fled object
## * select - allows one plot to be selected for printing. e.g. if you just want the plot for the fiber density set 'select=2'. When NULL (default) all  three plots are plotted on one page
## * density.scale - one of three options which define the scale on which the fiber/fine length densities should be plotted: "tree" (default) plots the densities of the fiber and fine length in a tree, "uncut.core" plots densities of uncut cells at least partially appeared in the increment core, "core" plots densities of observed (cut or uncut) cells in the increment core
{ if (is.null(ylab)) ylabel <- "density" else ylabel <- ylab
  
  ## plot functions for Kajaani data...
  plot.mixture <- function(x, rvec0, density.scale, xlab, main){
    if (is.null(xlab)) xlabel<- "cell length (mm)" else xlabel <- xlab
      if (is.null(main)){
         # if (x$model=="ggamma")  label<- "Generalized gamma mixture density"  
         # else     label<- "Log normal mixture density" 
         label<- "Fiber and fine density in the core"
      } else label <- main
      if (is.null(rvec0)) rvec <- seq(0,8-1e-3,length=300)
    ##  if (x$model=="lognorm") lx <- lx.logN
    ##  pred <- with(as.list(x$logpar),lx(x=rvec,theta=x$logpar,r=x$r,log=FALSE))
      if (density.scale=="tree")
        pred <- with(as.list(x$par), dw.mixture.withMu( x=rvec, par=x$par, r=x$r, mu.fibers=x$mu.fibers, mu.fines=x$mu.fines, prop.fines=x$prop.fines, model=x$model))
      else if (density.scale=="uncut.core")
        pred <- with(as.list(x$par), dy.mixture(x=rvec, par=x$par, model=x$model))
      else if (density.scale=="core")
         pred <- with(as.list(x$par), dx.mixture(x=rvec,par=x$par,r=x$r, model=x$model)) 
      hist(x$data,breaks=50,col="gray",freq=FALSE, xlab=xlabel, ylab=ylabel, ylim=c(0,max(pred)), main=label)      
      lines(rvec,pred,col=col,lwd=lwd)
  }

  plot.pdf.fiber <- function(x, rvec0, density.scale, xlab, main){
     if (is.null(xlab)) xlabel<- "fiber length (mm)" else xlabel <- xlab
     if (is.null(main)) {
         if (density.scale=="uncut.core")
             label <- "Density of fibers that at least partially appear in the core"
         else label<- substitute(paste("Fiber density in the ", density.scale),list(density.scale=density.scale)) 
     } else label <- main
     if (is.null(rvec0)) rvec <- seq(1e-14,8-1e-5,length=300)
     if (x$model=="ggamma")  par0 <- x$par[5:7] else if (x$model=="lognorm") par0 <- x$par[4:5]                   
     if (density.scale=="tree")
        pred <- with(as.list(par0), dw.fibers.withMu(x=rvec, par=par0, r=x$r, mu.fibers=x$mu.fibers, model=x$model))
     else if (density.scale=="uncut.core")
        pred <- with(as.list(par0), dy.fibers(x=rvec, par=par0, model=x$model))
     else if (density.scale=="core")
         pred <- with(as.list(par0), dx.fibers(x=rvec, par=par0, r=x$r, model=x$model)) 
     plot(rvec,pred,type="l", xlab=xlabel, ylab=ylabel, main=label, col=col+2, lwd=lwd)
  }

  plot.pdf.fine <- function(x, rvec0, density.scale, xlab, main){
     if (is.null(xlab)) xlabel<- "fine length (mm)" else xlabel <- xlab
     if (is.null(main)) label<- substitute(paste("Fine density in the ", density.scale),list(density.scale=density.scale))   else label <- main
     if (x$model=="ggamma"){
          if (is.null(rvec0)) rvec <- seq(0,2.5,length=300)    
          par0 <- x$par[2:4]         
      } else if (x$model=="lognorm"){ 
               if (is.null(rvec0)) rvec <- seq(1e-8,2.5,length=300) 
               par0 <- x$par[2:3]             
            }
     if (density.scale=="tree")
        pred <- with(as.list(par0), dw.fibers.withMu(x=rvec, par=par0, r=x$r, mu.fibers=x$mu.fines, model=x$model))
     else if (density.scale=="uncut.core")
        pred <- with(as.list(par0), dy.fibers(x=rvec, par=par0, model=x$model))
     else if (density.scale=="core")
         pred <- with(as.list(par0), dx.fibers(x=rvec, par=par0, r=x$r, model=x$model) )
     plot(rvec,pred,type="l", xlab=xlabel, ylab=ylabel, main=label, col=col+2, lwd=lwd)
  }

  ## plot functions for microscopy data...
  plot.microscopy.observed <- function(x, rvec0, xlab, main){
    if (is.null(xlab)) xlabel<- "fiber length (mm)" else xlabel <- xlab
      if (is.null(main)){
         label<- "Density of uncut fibers in the core"
      } else label <- main
      if (is.null(rvec0)) rvec <- seq(1e-8,8-1e-5,length=300)
     # if (x$model=="lognorm") fmic.1 <- f.logN.mic.1
     # pred <- with(as.list(x$logpar),fmic.1(x=rvec,theta=x$logpar,r=x$r,log=FALSE))
      pred <- with(as.list(x$par), dx.fibers.micro(x=rvec,par=x$par,r=x$r, model=x$model))
      hist(x$data,breaks=30,col="gray",freq=FALSE, xlab=xlabel, ylab=ylabel, main=label,cex.main=1.1) ## ylim=c(0,max(pred)),      
      lines(rvec,pred,col=col,lwd=lwd)
  }

  plot.pdf.fiber.micro <- function(x, rvec0, density.scale, xlab, main){
     if (is.null(xlab)) xlabel<- "fiber length (mm)" else xlabel <- xlab
     if (is.null(main)) label<- substitute(paste("Fiber length density in the ", density.scale),list(density.scale=density.scale))  else label <- main
     if (is.null(rvec0)) rvec <- seq(1e-14,8-1e-5,length=300)
     par0 <- x$par                   
     if (density.scale=="tree")
        pred <- with(as.list(par0), dw.fibers.micro.withMu(x=rvec, par=par0, r=x$r, mu.fibers=x$mu.fibers, model=x$model))
     else if (density.scale=="uncut.core")
        pred <- with(as.list(par0), dy.fibers.micro(x=rvec, par=par0, model=x$model))
     else if (density.scale=="core")
         pred <- with(as.list(par0), dx.fibers.micro(x=rvec, par=par0, r=x$r, model=x$model)) 
     plot(rvec,pred,type="l", xlab=xlabel, ylab=ylabel, main=label, col=col+2, lwd=lwd)
  }

  ##-----------------------------
  rvec0 <- rvec
if (x$data.type=="ofa"){
  if (is.null(select)) { ## plotting all three 
      old.par<-par(mfrow=c(2,2))
      ## plotting mixture density...
      plot.mixture(x=x, rvec0=rvec0, density.scale="core", xlab=xlab, main=main)
      ## plotting fiber density in a standing tree (on W scale)...
      plot.pdf.fiber(x=x, rvec0=rvec0, density.scale=density.scale, xlab=xlab, main=main)
      ## plotting fine density...
      plot.pdf.fine(x=x, rvec0=rvec0, density.scale=density.scale, xlab=xlab, main=main) 
  } else if (select==1) { ## plotting only mixure density...
      old.par<-par(mfrow=c(1,1))
      ## plotting mixture density...
      plot.mixture(x=x, rvec0=rvec0, density.scale="core", xlab=xlab, main=main)
    } else if (select==2) { ## plotting only fiber density...
        old.par<-par(mfrow=c(1,1))
        ## plotting fiber density...
        plot.pdf.fiber(x=x, rvec0=rvec0, density.scale=density.scale, xlab=xlab, main=main)
      } else if (select==3) { ## plotting only fine density...
          old.par<-par(mfrow=c(1,1))
          ## plotting fine density...
          plot.pdf.fine(x=x, rvec0=rvec0, density.scale=density.scale, xlab=xlab, main=main)    
        } else stop("not sure what to plot. 'select' must be either 1, 2, 3 or NULL.")
} else if (x$data.type=="microscopy"){
     if (is.null(select)) { ## plotting two plots... 
        old.par<-par(mfrow=c(1,2))
        ## plotting observed fiber density...
        plot.microscopy.observed(x=x, rvec0=rvec0, xlab=xlab, main=main)
      #  plot.mixture(x=x, rvec0=rvec0, density.scale="core", xlab=xlab, main=main)
        ## plotting fiber density in a standing tree (on W scale)...
        plot.pdf.fiber.micro(x=x, rvec0=rvec0, density.scale=density.scale, xlab=xlab, main=main)
     } else if (select==1) { ## plotting only mixure density...
         old.par<-par(mfrow=c(1,1))
         ## plotting observed fiber density...
         plot.microscopy.observed(x=x, rvec0=rvec0, xlab=xlab, main=main)
     } else if (select==2) { ## plotting only fiber density...
         old.par<-par(mfrow=c(1,1))
         ## plotting fiber density...
         plot.pdf.fiber.micro(x=x, rvec0=rvec0, density.scale=density.scale, xlab=xlab, main=main)
      } else stop("not sure what to plot. 'select' must be either 1, 2, or NULL.")
  }

  par(old.par)
} ## end plot.fled




