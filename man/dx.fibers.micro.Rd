\name{dx.fibers.micro}
\alias{dx.fibers.micro}
\alias{dy.fibers.micro}
\alias{dw.fibers.micro}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Density functions of the fiber length based on microscopy data}
\description{ Functions to get values of the density functions of the fiber length on three different scales based on microscopy data, being the uncut fibres in the core: the uncut fibres in the increment core (\code{dx.fibers.micro}), as true fiber lengths of those fibres that at least partially appears in the increment core (\code{dy.fibers.micro}) and as the true fiber lengths in the standing tree (\code{dw.fibers.micro}).
}
\usage{
  dx.fibers.micro(x, par, r, model="ggamma")
  dy.fibers.micro(x, par, model="ggamma")
  dw.fibers.micro(x, par, r, model="ggamma")
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
   \item{x}{vector of fiber length}
   \item{par}{vector of parameters for fiber length distribution}
   \item{r}{radius of the increment core} 
   \item{model}{either \code{ggamma} (default) or \code{lognorm}}   
}
%\details{ 
%}
\value{
returns a vector of density values 

}

\author{
Sara \enc{Sjöstedt}{Sjostedt} de Luna, Konrad Abramowicz, Natalya Pya Arnqvist
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{dx.fibers}}, \code{\link{dy.mixture}}    
}

\examples{
 library(fiberLD)
 r <- 2.5
 par <- c(1.8,2.7,2.6)
 x <- seq(.01, 2*r-.01,length=100)
 f1 <- dy.fibers.micro(x, par)
 plot(x,f1,type="l")

 f2 <- dx.fibers.micro(x, par, r)

 ## getting the density in the tree that goes beyond the length 
 ## of the diameter, 2r,... 
 w <- seq(0,8,length=200)
 f3 <- dw.fibers.micro(w, par, r)
 plot(w,f3,type="l")

}










