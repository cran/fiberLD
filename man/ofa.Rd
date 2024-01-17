\name{ofa}
\alias{ofa}
\docType{data}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Example of increment core data}

\description{Simulated data of cell lengths of both cut and uncut fines and fibres from an increment core (as measured by means of an optical fiber analyzer). The data were simulated from a censored version of a mixture of the fine and fiber length distributions under the assumption that the true lengths of those cells (fines and fibers) that at least partially appear in the increment core follow generalized gamma distributions. Useful for illustrating the use of \code{fled()}. The values of the parameters of the mixture density are \code{par=c(0.32, 0.001, 0.2757,5.022, 2.4, 3.3, 1.5)}, and the radius of the core is \code{r=2.5mm}.
}
\usage{
data(ofa)
}
%- maybe also `usage' for other objects documented here.

\format{ \code{ofa} is a vector of 3000 cell lengths simulated from the density 

\code{dx.mixture(x, par, r=2.5, model="ggamma")} with parameters  

\code{par=c(0.32, 0.001, 0.2757, 5.022, 2.4, 3.3, 1.5)}. 
} 


%\details{The data frame \code{ofa} relates to ......
%}

%\source{
%The data are adapted from Svensson, \enc{Sjöstedt}{Sjostedt} de Luna and Bondesson (2006).
%}

%\examples{
%}
\references{
Svensson, I., \enc{Sjöstedt}{Sjostedt} de Luna, S., Bondesson, L. (2006). Estimation of wood fibre length distributions from censored data through an EM algorithm. \emph{Scandinavian Journal of Statistics}, \bold{33}(3), 503--522.
}








