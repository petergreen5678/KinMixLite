\name{size}
\alias{size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate and display total size of BN tables for a DNA mixture}
\description{Calculate and display total size of BN tables for a DNA mixture}
\usage{
size(mixture)
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns total size, typically to be printed by bespoke method}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(KinMixLite)
data(test2data)
data(NGMDyes)

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
size(mixD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
