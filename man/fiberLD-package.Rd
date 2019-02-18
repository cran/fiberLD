\name{fiberLD-package}
%\Rdversion{0.1-0}
\alias{fiberLD-package}
\docType{package}
\title{Fiber Length Determination
}
\description{
\code{fiberLD} provides functions for estimating tree fiber (tracheid) length distributions in the standing tree based on increment core samples. Two types of data can be used with the package,   increment core data measured by means of an optical fiber analyzer (OFA), e.g. such as the Kajaani Fiber Lab, or measured by microscopy. Increment core data analyzed by OFAs consist of the cell lengths of both cut and uncut fibres (tracheids) and fines (such as ray parenchyma cells) without being able to identify which cells are cut or if they are fines or fibres. The microscopy measured data consist of the observed lengths of the uncut fibres in the increment core.
A censored version of a mixture of the fine and fiber length distributions is proposed to fit the OFA data, under distributional assumptions. The package offers two choices for the assumptions of the underlying density functions of the true fiber (fine) lenghts of those fibers (fines) that at least partially appear in the increment core, being the generalized gamma and the log normal densities.

Maximum likelihood estimation is used for estimating the model parameters for both the OFA analyzed data and the microscopy measured data. In addition a stochastic version of the expectation-maximization method is provided to fit the log normal model to the increment core data analyzed by OFAs. Details about methods and data can be found in references. 
}
\details{
\tabular{ll}{
Package: \tab fiberLD\cr
Type: \tab Package\cr
License: \tab GPL (>= 2)\cr
}}
\author{
Sara \enc{Sjöstedt}{Sjostedt} de Luna, Konrad Abramowicz, Natalya Pya Arnqvist

Maintainer: Natalya Pya Arnqvist<nat.pya@gmail.com>

}
%\note{
% } 
\references{
Svensson, I., \enc{Sjöstedt}{Sjostedt} de Luna, S., Bondesson, L. (2006). Estimation of wood fibre length distributions from censored data through an EM algorithm. \emph{Scandinavian Journal of Statistics}, \bold{33}(3), 503--522.

Chen, Z. Q., Abramowicz, K., Raczkowski, R., Ganea, S., Wu, H. X., Lundqvist, S. O., \enc{Mörling}{Morling}, T., \enc{Sjöstedt}{Sjostedt} de Luna, S., Gil, M.R.G., Mellerowicz, E. J. (2016). Method for accurate fiber length determination from increment cores for large-scale population analyses in Norway spruce. \emph{Holzforschung}. Volume \bold{70}(9), 829--838.

}

