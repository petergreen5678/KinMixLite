\name{require.compiled}
\alias{require.compiled}
\title{Force compilation of all BNs in a DNA mixture}
\description{Scan all Bayes nets in mixture, and compile any that are not already compiled}
\usage{
require.compiled(mixture)
}
\arguments{
  \item{mixture}{A DNAmixture object}
}
\value{No value is returned, the function is called for its side effect}
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

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,dyes=list(NGMDyes),
	triangulate=FALSE,compile=FALSE)
replace.tables.for.UAF(mixD,40,compile=FALSE)
require.compiled(mixD)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
