\name{dx.fibers}
\alias{dx.fibers}
\alias{dy.fibers}
\alias{dw.fibers}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Density functions of the fiber lengths}
\description{ Functions to get values of the density functions of the fiber length on three different scales: as observed in the increment core corresponding to cut and uncut fiber lengths in the core (\code{dx.fibers}), as true fiber lengths of the fibres that at least partially appear in the increment core (\code{dy.fibers}) and as the true fiber lengths in the standing tree (\code{dw.fibers}).
}
\usage{
  dx.fibers(x, par, r, model="ggamma")
  dy.fibers(x, par, model="ggamma")
  dw.fibers(x, par, r, model="ggamma")
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
   \item{x}{vector of fiber lengths}
   \item{par}{vector of parameters for fiber length distribution}
   \item{r}{radius of the increment core} 
   \item{model}{either \code{ggamma} (default) or \code{lognorm}}   
}
%\details{ Prints out the model ..... 
%}
\value{
returns a vector of density values 

}

\author{
Sara \enc{Sjöstedt}{Sjostedt} de Luna, Konrad Abramowicz, Natalya Pya Arnqvist
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{dx.mixture}}, \code{\link{dy.mixture}}, \code{\link{dw.mixture}}    
}

\examples{
 library(fiberLD)
 r <- 2.5
 par <- c(1.8,2.7,2.6)
 x <- seq(.01, 2*r-.01,length=100)
 f1 <- dy.fibers(x, par)
 plot(x,f1,type="l")

 f2 <- dx.fibers(x, par, r)
 f3 <- dw.fibers(x, par, r)

 ## the same functions can be used for plotting fine densities...
 par.fines <- c(-2, .5)
 x <- seq(0, 2,length=100)
 f1.fines <- dy.fibers(x, par.fines, model="lognorm")
 plot(x,f1.fines,type="l")
 
 ## getting the distribution in the tree that goes beyond the length 
 ## of the diameter, 2r ... 
 w <- seq(0,8,length=200)
 f5 <- dw.fibers(w,par,r) ## fiber density in the tree
 par.fines <- c(.06,1.5,3.3)
 f2.fines <- dw.fibers(w, par.fines,r)
 par(mfrow=c(1,2))
 plot(w,f5,type="l", main="fiber density")
 plot(w,f2.fines,type="l", main="fine density")

}











