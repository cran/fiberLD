\name{plot.fled}
%\Rdversion{1.0}
\alias{plot.fled}
\title{FLED plotting
}
\description{The function
   takes a \code{fled} object produced by \code{fled()} and creates several density plots. When the data consists of cell lengths from the increment core measured by an optical fiber analyzer ("ofa"), the function creates a histogram of the given data values together with the estimated density of the mixture model and two separate plots of the estimated fiber and fine lengths densities in the standing tree. With a microscopy sample (consisting of the lengths of uncut fibers in the increment core) the function creates two plots, a histogram of the given data with the estimated density of lenghts of the uncut fibers in the increment core and the estimated fiber length density in the standing tree.}

\usage{
\method{plot}{fled}(x,select=NULL, density.scale="tree", rvec=NULL, xlab=NULL,
                    ylab=NULL,  main=NULL, col=4, lwd=2,  ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a \code{fled} object as produced by \code{fled()}.}

   \item{select}{ allows one plot to be selected for printing. e.g., if you just want the plot for the fiber length density set \code{select=2}. When 'NULL' (default) all  three plots are plotted. 
}
   \item{density.scale}{one of three options which define the scale on which the fiber/fine length densities should be plotted: "tree" (default) plots the estimated densities of the fiber/fine lengths in the tree, "uncut.core" plots densities of cell lengths of those cells that at least partially appear in the increment core, "core" plots densities of the observed (cut or uncut) cell lengths in the increment core } 
   \item{rvec}{ values of cell lengths used to get estimates of densities.}

  \item{xlab}{If supplied then this will be used as the x label for all plots.}

  \item{ylab}{If supplied then this will be used as the y label for all plots.}

  \item{main}{Used as title for plots if supplied.}

%  \item{ylim}{If supplied then this pair of numbers are used as the y limits for each plot.}

 % \item{xlim}{If supplied then this pair of numbers are used as the x limits for each plot.}

  \item{col}{defines the color used for density plotting.}
  \item{lwd}{defines the line width.}

  \item{...}{ other graphics parameters to pass on to plotting commands.}

}

% \details{ to be added...
%%  ~~ If necessary, more details than the description above ~~
% }
\value{The function generates plots.

}

\author{
Sara \enc{Sjöstedt}{Sjostedt} de Luna, Konrad Abramowicz, Natalya Pya Arnqvist
}
\references{
 Svensson, I., \enc{Sjöstedt}{Sjostedt} de Luna, S., Bondesson, L. (2006). Estimation of wood fibre length distributions from censored data through an EM algorithm. \emph{Scandinavian Journal of Statistics}, \bold{33}(3), 503--522.

Chen, Z. Q., Abramowicz, K., Raczkowski, R., Ganea, S., Wu, H. X., Lundqvist, S. O., \enc{Mörling}{Morling}, T., \enc{Sjöstedt}{Sjostedt} de Luna, S., Gil, M.R.G., Mellerowicz, E. J. (2016). Method for accurate fiber length determination from increment cores for large-scale population analyses in Norway spruce. \emph{Holzforschung}. Volume \bold{70}(9), 829--838.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{  \code{\link{fled}}}

\examples{
## see ?fled help files
}


