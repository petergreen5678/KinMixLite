\name{pedigreeIBD}
\alias{pedigreeIBD}
\title{Construct IBD pattern distribution}
\description{Construct IBD pattern distribution from a pedigree and a target list of individuals}
\usage{
pedigreeIBD<-function (x, target = c("c1", "c2", "c3"), 
    cond = TRUE, quiet=TRUE, verbose = FALSE) 
}
\arguments{
  \item{x}{A pedigree in \code{\link{pedtools}} format}
  \item{target}{Character vector, some or all of the individual identifiers in the pedigree \code{x}}
  \item{cond}{should IBD pattern be condensed?}
  \item{quiet}{should resulting IBD pattern distribution be printed?}
  \item{verbose}{should trace information be printed?}
}

\details{
This function computes the multi-person condensed coefficients of identity for an arbitrary set of individuals, in the sparse notation of the IBD pattern distribution of Green and Vigeland (2019).
}

\value{IBD pattern distribution as a list with components \code{pr} and \code{patt}}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pedtools}}, \code{\link{KinMix}}
}

\references{
\emph{Multi-person condensed coefficients of identity},
by Peter J. Green and Magnus Dehli Vigeland,
University of Bristol technical report, 2019.
}

\examples{
require(KinMix)
require(ribd)

id<-c('gf','gm','b1','b2','m','c')
fid<-c(0,0,'gf','gf',0,'b1')
mid<-c(0,0,'gm','gm',0,'m')
sex<-c(1,2,1,1,2,0)   
x<-ped(id,fid,mid,sex)

IBD<-pedigreeIBD(x,c('m','c','b1','b2'))

kappaIBD(x,c('m','c','b1','b2'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
