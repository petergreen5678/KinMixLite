\name{rGTs}
\alias{rGTs}
\title{Simulate random genotype profiles and DNA samples for related individuals}
\description{Simulate random genotype profiles and DNA samples for arbitrarily related individuals}
\usage{
rGTs(nreps,IBD,db,DNA,sex=rep(0,ncontr),nU=0)
}
\arguments{
  \item{nreps}{Integer, number of replicates}
  \item{IBD}{Specification of relationships, as in \code{convertIBD}}
  \item{db}{Data frame, database of alleles and their frequencies, for each marker; see Formats section in \code{\link{KinMixLite}}.}
  \item{DNA}{Integer vector, numbers of DNA cells for the respective individuals, can be NA}
  \item{sex}{Integer vector, sex of the respective contributors: 1=male, 2=female, 0=unspecified}
  \item{nU}{Integer, include also this number of unrelated individuals}
}
\value{Data frame with variables \code{Sim}, \code{Sample.name}, \code{Marker}, \code{Allele}, and \code{DNA}, suitable for input to \code{simExtraction}, etc. See package \code{pcrsim}.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\details{
Genotype profiles are generated randomly, using the allele frequency database \code{db}, under the relationships specified by the \code{IBD} argument. In accordance with the underlying biology, allele values for the AMEL marker (if this is one of the markers included) are not influenced by relationships with other individuals; however they are influenced by the sex of the individuals, where this is known. Information on sex can be specified by the optional argument \code{sex}: a male is given the profile X-Y, a female X-X, and an individual with unspecified sex X-X or X-Y with equal probabilities.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pcrsim-package}}
}
\examples{
require(KinMixLite)
data(test2data)
data(NGMDyes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
