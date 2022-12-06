\name{loop.rpt.IBD}
\alias{loop.rpt.IBD}
\title{Analysis of DNA mixtures with familial relatioships}
\description{
Analysis of DNA mixtures with familial relationships, by looping over traces, markers, and IBD patterns, to reduce total BN table size, at some price in execution time}
\usage{
loop.rpt.IBD(listdata, pars, IBD, typed.gts = NULL, inds = 1, 
    jtyped = ncol(IBD$patt)/2 - length(typed.gts) + seq_along(typed.gts), 
    jcontr = seq_along(inds), targets = NULL, contribs, 
    quiet=FALSE, verbose=FALSE, presence.only=FALSE, ...) 
}
\arguments{
  \item{listdata}{as in call to DNAmixture}
  \item{pars}{parameter structure, in \code{mixpar} format}
  \item{IBD}{multi-person coefficients of identity, in any of the formats accepted by \code{convertIBD}}
  \item{typed.gts, inds, jtyped, jcontr, targets, contribs, quiet
	}{as in call to rpt.IBD}
  \item{verbose}{should per-marker and overall log10LR's be reported?}
  \item{...}{other arguments to DNAmixture, particularly including \code{k}, \code{C}, \code{database}}
}

\value{The value of the overall \code{log10 LR}, and the contributions of individual markers in the form of a vector-valued attribute 'log10LR', are returned invisibly; individual marker/pattern values are also printed out.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(KinMix)
data(test2data)
data(NGMDyes)

C<-0.001

## Fit 3-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db)
pars3<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.6,U2=0.3,U3=0.1)))
baseline3<-logL(mixD)(pars3)
size(mixD)

## Fit 3-person mixture - in which U1 and U2 have a parent-child relationship

mixD<-DNAmixture(list(epg),k=3,C=rep(list(C),length(list(epg))),database=db,
	triangulate=FALSE,compile=FALSE)
delete.DQnodes(mixD)
rpt.IBD(mixD,kin=c(0,1,0),typed.gts=list(),inds=1:2,jtyped=NULL)
size(mixD)
log10LR<-(logL(mixD)(pars3)-baseline3)/log(10)
cat('log10 LR',log10LR,'\n')

## the same analysis by loop.rpt.IBD

listdata<-list(epg)
loop.rpt.IBD(listdata,pars3,kinship=c(0,1,0),
	k=3,C=rep(list(C),length(listdata)),database=db,
	typed.gts=list(),inds=1:2,jtyped=NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
