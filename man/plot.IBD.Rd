\name{plot.IBD}
\alias{plot.IBD}
\title{Plot IBD patterns and pattern distributions}
\description{Plot IBD patterns and pattern distributions}
\usage{
plot.IBD(pattern,labels=NULL,
	probs=NULL,order=NULL,colrs=c('black','red','blue'),
	digits=3,nr=ceiling(sqrt(np))
}
\arguments{
  \item{pattern}{A matrix whose rows are IBD patterns, or a list whose components are \code{patt}, such a matrix, together with \code{pr}, a vector of the corresponding probabilities}
\item{labels}{Vector of numerical or character labels for the patterns, if \code{NA}, labels are constructed from the patterns by catenation, if \code{NULL}, the labels are not displayed.}
\item{probs}{Vector of probabilities of the patterns, if not provided as a component of \code{pattern}; if \code{NULL}, the probabilities are not displayed.}
\item{order}{A character string, partially matched using \code{pmatch} to one of 'pattern', 'probs', or 'labels', requesting ordering diagram accordingly (in the case of \code{probs} in decreasing order, \strong{or} a numeric, complex, character or logical vector of length the number of patterns, requesting ordering by this variable, \strong{or} NULL (the default), requesting no re-ordering.}
\item{colrs}{A vector of colours: ties in the ordering variable are indicated by coloured groups, with colours chosen cyclically from this vector.}
\item{digits}{Integer, overwriting default number of significant digits for \code{probs}}
\item{nr}{Integer, overwriting default number of rows for plotted array, default a rounding up of the square root of the number of patterns.}
}
\value{No value is returned, the function is called for its side effect, a plot on the current display device.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(KinMixLite)
require(ribd)
data(emperors)

plot.IBD(convertIBD('3cousins-star'),order='probs',col=c('blue','red','black'))

plot(attr(emperors,'ped'))

o<-order(emperors$pr)[1:12]
plot.IBD(emperors$patt[o,],probs=emperors$pr[o],labels=NA,order='probs')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
