\name{add.child.meiosis.nodes}
\alias{add.child.meiosis.nodes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Replace CPTs for mixture contributor a putative father, given Child genotype}
\description{loop over markers, and alleles within markers
to create nodes for child allele count nodes,
for paternity model with only Child genotyped
then compile all domains. Implements method MBN.
}
\usage{
add.child.meiosis.nodes(mixture,aca,ind=1)
}
\arguments{
  \item{mixture}{A compiled DNAmixture object}
  \item{aca}{Child's genotype profile as an allele count array}
  \item{ind}{Index of contributor regarded as Parent (or Child): which `unknown' contributor are we modelling by amending his/her CPTs?}
}
\details{
To calculate the likelihood of this model, conditional on the child's genotype, a call to this function should be followed by (a) setting the finding of the child's genotype by defining \code{extra.findings}, (b) evaluating the loglikelihood using \code{logLX}, and (c) correcting the result by subtracting the log probability of the child's genotype, all as in the example below. Without (c), the value returned is that of the peak heights \emph{and} the child's genotype.
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

mixMBN<-DNAmixture(list(epg),k=2,C=list(C),database=db,
	dyes=list(NGMDyes),triangulate=F,compile=F)
cgtcaca<<-gt2aca(mixMBN,Cgt)
add.child.meiosis.nodes(mixMBN,cgtcaca,1)
extra.findings<-make.findings(list(
	list('Male',ind=1),
	list('Caca',aca='cgtcaca')
	))
likelihood<-exp(logLX(mixMBN)(pars)-attr(cgtcaca,'logGt'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


