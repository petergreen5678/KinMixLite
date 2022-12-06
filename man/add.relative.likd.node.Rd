\name{add.relative.likd.node}
\alias{add.relative.likd.node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replace CPTs for mixture contributor a putative father, given Child genotype}
\description{loop over markers, and alleles within markers
to create node Rlikd for relative likelihood for individual i,
for paternity model with only Child genotyped
then compile all domains. Implements method ALN.
}
\usage{
add.relative.likd.node(mixture,aca,ind=1)
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
  \item{aca}{Child's genotype profile as an allele count array}
  \item{ind}{Index of contributor regarded as Parent (or Child): which `unknown' contributor are we modelling by amending his/her CPTs?}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
require(KinMix)
data(test2data)
data(NGMDyes)

# set threshold C
C<-0.001

pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.9,U2=0.1)))

mixALN<-DNAmixture(list(epg),k=2,C=list(C),database=db,
	dyes=list(NGMDyes),triangulate=F,compile=F)
cgtcaca<<-gt2aca(mixALN,Cgt) 
add.relative.likd.node(mixALN,cgtcaca,1)
extra.findings<-make.findings(list(
	list('Male',ind=1),
	list('Rlikd',aca='cgtcaca',cgt='Cgt',evid='Revid')
	))
likelihood<-exp(logLX(mixALN)(pars))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
