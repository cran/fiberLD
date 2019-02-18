\name{dx.mixture}
\alias{dx.mixture}
\alias{dy.mixture}
\alias{dw.mixture}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Mixture density functions of the cell lengths in the increment core}
\description{ Functions to get values of the mixture density functions of the cell lengths on three different scales: as observed in the increment core, i.e. cut or uncut fibers or fines (\code{dx.mixture}), as true cell lengths (fibers or fines) that at least partially appears in the increment core (\code{dy.mixture}) and as the true cell lengths (fines or fibres) in the standing tree (\code{dw.mixture}).
}
\usage{
  dx.mixture(x, par, r, model="ggamma")
  dy.mixture(x, par, model="ggamma")
  dw.mixture(x, par, r, model="ggamma")
}
%- maybe also `usage' for other objects documented here.
\arguments{ 
   \item{x}{vector of cell length values}
   \item{par}{vector of mixture model parameters}
   \item{r}{radius of the increment core} 
   \item{model}{either \code{ggamma} (default) or \code{lognorm}}   
}

\value{
returns a vector of density values 

}

\author{
Sara \enc{Sjöstedt}{Sjostedt} de Luna, Konrad Abramowicz, Natalya Pya Arnqvist
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{dx.fibers}}, \code{\link{dy.fibers}}, \code{\link{dw.fibers}}    
}

\examples{
\dontrun{ 
 library(fiberLD)
 data(cell.length) 
 r <- 6
 d <- fled(data=cell.length,model="lognorm",r=r)
 x <- seq(0, 8,length=100)
 f1 <- dx.mixture(x=x, par=d$par,r=r, model="lognorm")
 plot(x,f1,type="l")

 w <- seq(0, 2*r+2,length=150)
 f2 <- dy.mixture(w, d$par, model="lognorm")
 f3 <- dw.mixture(w, d$par, 6, model="lognorm")
 plot(w,f3,type="l")
}
}











