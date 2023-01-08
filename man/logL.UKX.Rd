\name{logL.UKX}
\alias{logL.UKX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Evaluates mixture log likelihood for unknown contributors}
\description{Replacement for logL.UK in DNAmixtures that calls \code{extra.findings} immediately before 
propagating all findings and returning the normalising constant for the network.}
\usage{logL.UKX(mixture, initialize = TRUE)}
\arguments{
  \item{mixture}{Compiled DNAmixture object.}
  \item{...}{All as in \code{mixML}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{The log likelihood.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{logL.UK}}.
}
\examples{
require(KinMixLite)
data(test2data)
data(NGMDyes)

# set threshold C
C<-0.001

pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.9,U2=0.1)))

mixMBN<-DNAmixture(list(epg),k=2,C=list(C),database=db,
	dyes=list(NGMDyes),triangulate=F,compile=F)
cgtcaca<<-gt2aca(mixMBN,Cgt)
add.child.meiosis.nodes(mixMBN,cgtcaca,1)
extra.findings<-make.findings(list(
	list('Male',ind=1),
	list('Caca',aca='cgtcaca')
	))

logL.UKX(mixMBN)(pars)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
