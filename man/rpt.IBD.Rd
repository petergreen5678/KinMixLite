\name{rpt.IBD}
\alias{rpt.IBD}
\alias{rpt.typed.relatives}
\alias{rpt.typed.child}
\alias{replace.Ui.tables}
\alias{rpt.typed.parents}
\alias{rpt.typed.relative}
\title{Replacing CPTs for selected mixture contributors with familial relationships}
\description{
Used after a call to \code{DNAmixture} with \code{compile=FALSE,triangulate=FALSE}, this function replaces the CPTs for the genotype allele count arrays for specified mixture contributors by those representing the specified relationship with each other and typed relatives
}
\usage{
rpt.IBD(mixture, IBD, typed.gts = NULL, inds = 1, 
jtyped = ncol(IBD$patt)/2 - length(typed.gts) + seq_along(typed.gts), 
jcontr = seq_along(inds), 
targets=NULL, contribs, quiet=FALSE, all.freq = NULL, compile = TRUE) 

rpt.typed.relatives(mixture, IBD, typed.gts = NULL, inds = 1, 
jtyped = ncol(IBD$patt)/2 - length(typed.gts) + seq_along(typed.gts), 
jcontr = seq_along(inds), 
targets=NULL, contribs, quiet=FALSE, all.freq = NULL, compile = TRUE)

rpt.typed.child(mixture, Cgt, ind=1)

replace.Ui.tables(mixture, Cgt, ind=1)

rpt.typed.parents(mixture, Mgt, Fgt, ind=1, compile=TRUE)

rpt.typed.relative(mixture, Rgt, IBD=c(0.25,0.5,0.25), ind=1, compile=TRUE)
}
\arguments{
  \item{mixture}{DNAmixtures object created by previous call to \code{DNAmixture} with \cr
	\code{triangulate=FALSE,compile=FALSE}}
  \item{IBD}{relationships between the specified individuals, as multi-person condensed 
coefficients of identity, in one of several representation; see Details.}
  \item{typed.gts}{list of 0 or more genotypes of relatives; the components of this list must be named (with the id's of the relevant individuals) if \code{targets} and \code{contribs} are used to code the correspondences (see Details).}
  \item{inds}{vector of 1 or more integers: which `unknown' contributors are we modelling by amending their CPTs? The elements should be listed in the same order as the corresponding pairs of columns of the IBD patterns in \code{IBD}}
  \item{jtyped}{indices of pairs of columns of \code{IBD$patt} that correspond to the 
typed relatives (if any); default the last \code{length(typed.gts)} pairs of columns}
  \item{jcontr}{indices of pairs of columns of \code{IBD$patt} that correspond to the 
relevant mixture contributors; default the first \code{length(inds)} pairs of columns}
  \item{targets}{Character vector of the tags of the individuals referred to in \code{IBD}}
  \item{contribs}{Character vector of the tags of the individuals included in the mixture,
in order}
  \item{quiet}{should calculated values of inds, jtyped and jcontr be reported?}
  \item{all.freq}{alternative allele frequency database(s), see Details.}
  \item{compile}{logical flag: should mixture object be compiled on exit?}
  \item{ind}{as \code{inds}, used above when only one allowed}
  \item{Cgt, Mgt, Fgt, Rgt}{individual genotypes}
}
\details{
In using \code{rpt.IBD} or \code{rpt.typed.relatives} (which is identical), the correspondence between mixture contributors, specified relationships, and typed genotype profiles should be specified \strong{either} with \code{inds}, \code{jcontr} and \code{jtyped}, \strong{or} using \code{targets}, \code{contribs} and through the names of the components in \code{typed.gts}: the two representations should not be mixed up.

Special cases are treated slightly more efficiently:
\code{rpt.typed.child}: single contributor, single typed relative, parent or child;
\code{rpt.typed.parents}: single contributor, both parents typed;
\code{rpt.typed.relative}: single contributor, single typed relative.

Note that \code{IBD$patt} always has an even number of columns, two for each individual in the joint relationship specified; \code{jtyped} and \code{jcontr} are vectors of indices of these individuals, i.e. to pairs of adjacent columns of \code{IBD$patt}.

Multiple functions in this group can validly be called sequentially (with all but the last having \code{compile=FALSE}) providing they reference different sets of contributors.

There are multiple valid representations for relationships in the argument \code{IBD} -- as an IBD pattern distribution, via a pedigree, or. in the case of just two individuals. via either a vector of 3 kappas or 9 Deltas (Jacquard's condensed coefficients of identity). For full details, see \code{\link{convertIBD}}.
If \code{IBD} is missing, the default value represents parent-child.

In the interests of upward compatibility, in \code{rpt.typed.child} and \code{replace.Ui.tables} (which are identical), the argument \code{Cgt} can be given as either a genotype profile data frame, or an allele count array.

By default, the allele frequency database used for the founding genes is that used when the \code{mixture} object is created, in an earlier call to \code{DNAmixture}. A non-null value for the \code{all.freq} argument allows the user to specify alternative database(s) for the founding genes. If its value is an allele frequency database (in the format specified in \code{\link{KinMix}}) then that database is used for all founding genes; if the value of the argument is a list of such databases, then component k of the list is used for allele frequencies for the founding gene labelled k in the \code{IBD} argument. Note that this option allows modelling of mixtures where different contributors are drawn from different populations, whether or not there are relationships among individuals. 
}

\value{Vector of marker-specific probabilities of the typed genotypes.}
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

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

## Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,,list(Cgt)) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture model in which contributor 1 is father of a typed individual Cgt 
## with mother Mgt

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,,list(Mgt,Cgt)) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 2-person mixture model in which contributors are two parents of a child with 
## genotype Cgt, and a parent of one of them has genotype Rgt. Note the encoding of allele 
## labels to reduce the complexity of the IBD pattern distribution IBD.

IBD<-list(patt=rbind(c(1,3,2,4,1,2,1,5),c(1,3,2,4,1,2,3,5)))

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,IBD,list(Cgt,Rgt),1:2) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## the same, with individuals and relationships denoted by character tags

mixD<-DNAmixture(list(epg),k=2,C=list(0.001),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,IBD,list(c=Cgt,gf=Rgt),targets=c('f','m','c','gf'),contribs=c('f','m')) 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
