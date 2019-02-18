\name{fled}
%\Rdversion{1.0}
\alias{fled}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fiber length determination
}
\description{This function estimates fiber (tracheid) and fine (e.g. ray parenchyma cells and other small particles) lengths distribution in standing trees based on increment cores (cylindric wood samples). The data from the increment cores contain uncut fiber, fibers cut once or twice (cut by the borer) as well as non-fiber cells so-called 'fines'. A censored version of a mixture of the fine and fiber length distributions is therefore proposed to fit the data. The function offers two choices for the underlying density functions of the true unobserved uncut lengths of the fines and fibers in the increment core such as generalized gamma and log normal densities. The parameters of the mixture models are estimated by log likelihood maximization. The routine calls an \code{optim()} or \code{nlm()} functions for optimization procedure with the possibility to use a supplied gradient function. Some parameters of the generalized gamma mixture model can be fixed (rather than estimated) at the given values. 
}
\usage{
 fled(data=stop("No data supplied"), data.type="ofa", r=2.5, 
     model="ggamma", method="ML", parStart=NULL, fixed=NULL, 
    optimizer=c("optim","L-BFGS-B","grad"),lower=-Inf,upper=Inf,...)
}
\arguments{
   \item{data}{A numeric vector of cell lengths from increment cores.}
   \item{data.type}{type of data supplied: ''ofa'' (default) measured by an optical fiber analyser, or measured by ''microscopy'' (only the lengths of uncut fibers in the core). }
   \item{r}{radius of the increment core (default 2.5).} 
   \item{model}{if \code{model="ggamma"} then the distributions of the true lengths of the fibers (fines) that at least partially appear in the increment core are assumed to follow generalized gamma distributions; if \code{model="lognorm"} then log normal distributions are assumed on those fiber (fine) lengths.}
  \item{method}{Currently only maximum likelihood method' 'ML'' is available.}
  \item{parStart}{numerical vector of starting values of parameters (or fixed values for ggamma model when \code{!is.null(fixed)}). The parameter values of the generalized gamma model should be given in the following order, 
 
\eqn{(\epsilon, b_{fines},d_{fines},k_{fines},b_{fibers},d_{fibers},k_{fibers})}.

The parameter values of the log normal model are in the order   

\eqn{(\epsilon, \mu_{fines}, \sigma_{fines}, \mu_{fibers}, \sigma_{fibers})} (see Details below). }
   \item{fixed}{ TRUE/FALSE vector of seven components used to tell which parameters of ggamma model to fix. These are fixed at the values given in the argument \code{parStart}).    The positive values in \code{parStart} for non-fixed parameters are starting values for the optimiser, the negative or zero values indicate that no starting values are assumed. Note, fixing parameter values currently works only with 'optim'.}
   \item{optimizer}{numerical optimization method used to minimize 'minus' the loglikelihood function of the observed data: 'optim', 'nlm' or 'nlm.fd' (nlm is based on finite-difference approximation of the derivatives). If optimizer==''optim'' then the second argument specifies the numerical method to be used in 'optim' (''Nelder-Mead'', ''BFGS'', ''CG'', ''L-BFGS-B'', ''SANN''. The third element of 
\code{optimizer} indicates whether the finite difference approximation should be used ('fd') or analytical gradient ('grad') for the 'BFGS', 'CG' and 'L-BFGS-B' methods. The default is \code{optimizer=c("optim",} \code{"L-BFGS-B","grad")}. } 
   \item{lower, upper}{Bounds on the parameters for the "L-BFGS-B" method. The order of the bounds values has to be the same as the order of the \code{parStart}. Note that these bounds are on the original rather than transformed scale of the parameters used for optimization.} 
   \item{...}{Further arguments to be passed to \code{optim}. }
}
\details{
  The probability density function of the three-parameter generalized gamma distribution proposed by Stacy (1962) can be written as
  \deqn{f(y;b,d,k) = d   b^{-d k}   y^{d k-1}  \exp[-(y/b)^d] / \Gamma(k),}{f(y;b,d,k) = d * b^(-d*k) * y^(d*k-1) * exp(-(y/b)^d)/gamma(k),}
  where \eqn{b > 0}, \eqn{d > 0}, \eqn{k > 0},
  and \eqn{y > 0}.

The probability density function of the log normal distribution can be written as
  \deqn{f(y;\mu, \sigma) =\exp[-(\log (y)-\mu)^2/(2\sigma^2)]/(y \sigma\sqrt{2\pi}), }{f(y;mu, sigma) = exp(-(log(y)-mu)^2/(2*sigma^2))/(y * sigma * sqrt(2*pi)),}
 where \eqn{\sigma > 0} and \eqn{y > 0}. 

}


\value{
  \item{cov.par}{approximate covariance matrix of the estimated parameters.}
  \item{cov.logpar}{approximate covariance matrix of the transformed estimated parameters.}
  \item{loglik}{the log likelihood value corresponding to the estimated parameters.}
  \item{model}{model used}
  \item{mu.fibers}{estimated mean value of the fiber lengths in the standing tree.}
  \item{mu.fines}{estimated mean value of the fine lengths in the standing tree.}
  \item{mu.cell}{estimated mean value of the cell lengths in the standing tree.}
  \item{prop.fines}{estimated proportion of fines in the standing tree.}
  \item{par}{the estimated parameters on the original scale.}
  \item{logpar}{the estimated values of the transformed parameters.}
  \item{termcode}{an integer indicating why the optimization process terminated (see \code{\link{optim}}).}
  \item{conv}{indicates why the optimization algorithm terminated.}
  \item{iterations}{number of iterations of the optimization method taken to get convergence.}
  \item{fixed}{TRUE/FALSE vector denoting if a parameter of ggamma model is fixed or not.}
  \item{n}{number of observations}
}

\note{ 
  The idea and some of the code for fixing parameters with \code{optim()} is due to Barry Rowlingson, October 2011.

}

\author{
Sara \enc{Sjöstedt}{Sjostedt} de Luna, Konrad Abramowicz, Natalya Pya Arnqvist
}
\references{
Svensson, I., \enc{Sjöstedt}{Sjostedt} de Luna, S., Bondesson, L. (2006). Estimation of wood fibre length distributions from censored data through an EM algorithm. \emph{Scandinavian Journal of Statistics}, \bold{33}(3), 503--522.

Chen, Z. Q., Abramowicz, K., Raczkowski, R., Ganea, S., Wu, H. X., Lundqvist, S. O., \enc{Mörling}{Morling}, T., \enc{Sjöstedt}{Sjostedt} de Luna, S., Gil, M.R.G., Mellerowicz, E. J. (2016). Method for accurate fiber length determination from increment cores for large-scale population analyses in Norway spruce. \emph{Holzforschung}. Volume \bold{70}(9), 829--838.


Stacy, E. W. (1962). A generalization of the gamma distribution.
  \emph{Annals of Mathematical Statistics}, \bold{33}(3), 1187--1192.

}

\section{Warning }{Fixing the parameters with the generalized gamma model may lead to unstable results of the optim method.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%    \code{\link{...}}, \code{}    
%}

\examples{

library(fiberLD)
## using microscopy data (uncut fiber lengths in the increment core)
data(microscopy)
dat <- microscopy[1:200]
m1 <- fled(data=dat,data.type="microscopy",model="ggamma",r=2.5) 
summary(m1)
plot(m1)

## and with log normal model...
m2 <- fled(data=dat,data.type="microscopy",model="lognorm",r=2.5)
summary(m2)
plot(m2)

\dontrun{ 
## using data measured by an optical fiber analyser
data(cell.length) 
d1 <- fled(data=cell.length,model="lognorm",r=6)
summary(d1)
plot(d1)
x11()
plot(d1,density.scale="uncut.core")

## change the model to generalized gamma
## and set lower and upper bounds on the parameters for 
## the "L-BFGS-B" method ... 
d2 <- fled(data=cell.length,model="ggamma",r=6,lower=c(.12,1e-3,.05,rep(.3,4)),
      upper=c(.5,2,rep(7,5)))  
d2
summary(d2)
plot(d2,select=1)
}

}


