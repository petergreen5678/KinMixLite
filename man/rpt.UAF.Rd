\name{rpt.UAF}
\alias{rpt.UAF}
\alias{replace.tables.for.UAF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{replace.tables.for.UAF
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replace.tables.for.UAF(mixture, M, compile = TRUE)

rpt.UAF(mixture, M, compile = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mixture}{DNAmixtures object created by previous call to \code{DNAmixture} with \cr
\code{triangulate=FALSE,compile=FALSE}}
  \item{M}{Size of allele frequency database}
  \item{compile}{logical flag: should mixture object be compiled on exit?}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{No value is returned, the function is called for its side effect}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\note{
%%  ~~further notes~~
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
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

## Fit 2-person mixture model under assumption that database size was only 40

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
replace.tables.for.UAF(mixD,40)
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
