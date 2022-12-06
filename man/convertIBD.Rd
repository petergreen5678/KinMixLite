\name{convertIBD}
\alias{convertIBD}
\title{Construct IBD pattern distribution}
\description{Construct IBD pattern distribution from one of several alternative representations
of multi-person condensed coefficients of identity}
\usage{
convertIBD(x)
}
\arguments{
  \item{x}{A string, a vector of length 3 or 9, a list with components \code{pr} and \code{patt}, or a list with two components, a pedigree and a vector of target id's; see Details}
}

\details{
Possible formats for the input \code{x} are:
\enumerate{
\item certain verbal mnemonics; currently one of the following (or an unambiguous partial match):
c('sibs','parent-child','half-sibs',
	'cousins','half-cousins','second-cousins',
	'double-first-cousins',
	'quadruple-half-first-cousins',
	'3cousins-cyclic','3cousins-star','trio')
\item a vector of 3 kappas
\item a vector of 9 Deltas
\item a list with matrix or vector valued component \code{patt}, with or without component \code{pr}
\item a list with 2 components, the first being a pedigree in the sense of the \code{pedtools} package, the second a vector of target id's
}
}

\value{IBD pattern distribution as a list with components \code{pr} and \code{patt}}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pedtools}}, \code{\link{KinMix}}


}
\examples{
require(KinMix)
data(test2data)
data(NGMDyes)

IBD<-convertIBD('parent-child')

IBD<-convertIBD(c(0.5,0.5,0.0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
