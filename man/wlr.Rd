\name{wlr}
\alias{wlr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes paternity LR using WLR method}
\description{Computes overall LR from Ugt-specific LR's using estimated Ugt genotype
profile in \code{sep} corresponding to contributor i in the mixture as Father; uses 
Child genotype information in \code{Cgt} data.frame and optionally Mother's genotype in \code{Mgt}.
Implements method WLR.
}
\usage{
wlr<-function(sep, Cgt, db, ind=1, Mgt=NULL)
}
\arguments{
  \item{sep}{Separation, a list of configurations of genotypes for some or all unknown contributors, output by \code{map.genotypes}.}
  \item{Cgt}{Child's genotype profile as a data frame containing variables \code{marker}, \code{allele1} and \code{allele2}}
  \item{db}{Allele frequency database}
  \item{ind}{Index of contributor regarded as Father}
  \item{Mgt}{(optionally) Mother's genotype profile as a data frame containing variables \code{marker}, \code{allele1} and \code{allele2}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns LR for paternity}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{map.genotypes}}.
}
\examples{
require(KinMixLite)
data(test2data)
data(NGMDyes)

# set threshold C
C<-0.001

pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.9,U2=0.1)))

mixWLR<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes))
setPeakInfo(mixWLR,pars)
sepWLR<-map.genotypes(mixWLR,type="all",pmin=0.0001,U=1)
LR<-wlr(sepWLR,Cgt,db)
cat('\nWLR LR:',LR,'; log10(LR):',log10(LR),'\n')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
