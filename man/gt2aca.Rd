\name{gt2aca}
\alias{gt2aca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converts genotype profile format}
\description{Returns list of vectors of allele counts corresponding to genotype profile in gt}
\usage{
gt2aca<-function(mixture,gt,eps=0)
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
  \item{gt}{Genotype profile as a data frame containing variables \code{marker}, \code{allele1} and \code{allele2}}
  \item{eps}{If non-zero, the function creates the output allele count arrays in a different format, that mitigates subsequent propagation errors in some situations. Instead of a vector of allele counts, each element of the list is a matrix with 3 columns, corresponding to allele counts 0, 1 and 2, with entries 1 or \code{eps}. This format is equally acceptable for the \code{aca} parameter of the \code{Aca} option in the \code{\link{make.findings}} function.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns list of vectors of allele counts. The log probability for the genotype is returned in its attribute '\code{logGt}'.}
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

# set threshold C
C<-0.001

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes))
cgtcaca<<-gt2aca(mixD,Cgt)

print(Cgt)
print(cgtcaca)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
