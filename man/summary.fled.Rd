\name{summary.fled}
%\Rdversion{1.1-0}
\alias{summary.fled}
\alias{print.summary.fled}

\title{Summary for a fled fit
}
\description{Takes a \code{fled} object produced by \code{fled()} and produces various useful
summaries from it. 
}
\usage{
\method{summary}{fled}(object,...)

\method{print}{summary.fled}(x,digits = max(3, getOption("digits") - 3),...)
}

\arguments{
\item{object}{a fitted \code{fled} object as produced by \code{fled()}.}

\item{x}{a \code{summary.fled} object produced by \code{summary.fled()}.} 

\item{digits}{controls the number of digits printed in the output.}

\item{...}{ other arguments.}

}
% \details{
%%  ~~ If necessary, more details than the description above ~~
% }Filter


\value{
\code{summary.fled} produces the following list of summary information for a \code{fled} object. 


\item{fixed}{is a logical vector for any model parameters that are fixed}
\item{n}{number of observations}
\item{loglik}{minimized minus log likelihood for the model}
\item{model}{the model used}
\item{method}{'ML' method used}
\item{data.type}{type of data used}
\item{conv}{indicates why the optimization algorithm terminated}
\item{p.table}{table of model parameters}
\item{ss.table}{table of summary statistics for cell lengths in the increment core}
\item{w.fine}{table of summary statistics for fine lengths in the standing tree}
\item{w.fiber}{table of summary statistics for fiber lengths in the standing tree}
\item{mean.w}{expected value of the cell lengths in the standing tree}
\item{eps.tree}{proportion of fines in the standing tree}
\item{se.eps.tree}{standard error of eps.tree}
}


\author{
Sara \enc{Sjöstedt}{Sjostedt} de Luna, Konrad Abramowicz, Natalya Pya Arnqvist
}
\references{
Svensson, I., \enc{Sjöstedt}{Sjostedt} de Luna, S., Bondesson, L. (2006). Estimation of wood fibre length distributions from censored data through an EM algorithm. \emph{Scandinavian Journal of Statistics}, \bold{33}(3), 503--522.

Chen, Z. Q., Abramowicz, K., Raczkowski, R., Ganea, S., Wu, H. X., Lundqvist, S. O., \enc{Mörling}{Morling}, T., \enc{Sjöstedt}{Sjostedt} de Luna, S., Gil, M.R.G., Mellerowicz, E. J. (2016). Method for accurate fiber length determination from increment cores for large-scale population analyses in Norway spruce. \emph{Holzforschung}. Volume \bold{70}(9), 829--838.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{fled}} 
}

\examples{
## see ?fled help files
}

\keyword{fiber length} \keyword{mixture model} \keyword{generalized gamma}


