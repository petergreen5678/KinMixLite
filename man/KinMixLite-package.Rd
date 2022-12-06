\name{KinMixLite-package}
\alias{KinMixLite-package}
\alias{KinMixLite}
\docType{package}
\title{
\packageTitle{KinMixLite}
}
\description{
\packageDescription{KinMixLite}
}
\details{
This package is a toolkit for inference about mixtures and familial relationships, either between contributors or between a contributor and other typed individuals. It extends the functionality proposed in Green and Mortera (2017) by allowing more general relationships, specified in general by an IBD pattern distribution - the generalisation to more than two individuals of the coefficients of identity of Jacquard (1974). Details are in a paper in preparation by Green and Mortera. \code{KinMixLite} extends the capability of the \pkg{DNAmixturesLite} package, and intimately relies on that package; it again uses the \pkg{RHugin} package for Bayes Net calculations.
}
\author{
\packageAuthor{KinMixLite}

Maintainer: \packageMaintainer{KinMixLite}
}
\references{
Green, P. J. and Mortera, J. (2017). Paternity testing and other inference about relationships from DNA mixtures.
\emph{Forensic Science International: Genetics}. \url{http://dx.doi.org/10.1016/j.fsigen.2017.02.001}. 

Jacquard, A. (1974) 
\emph{The genetic structure of populations}. Springer-Verlag.
}
\keyword{ package }
\seealso{
\code{\link[DNAmixturesLite:DNAmixturesLite-package]{DNAmixturesLite}}
}

\section{Formats}{
An \bold{allele frequency database} is a data frame containing variables \code{marker}, \code{allele} and \code{frequency} (character, numeric and numeric respectively).

A \bold{mixture profile} is a data frame containing variables \code{marker}, \code{allele} and \code{height} (character, numeric and numeric respectively).

A \bold{genotype profile} is a data frame containing variables \code{marker}, \code{allele1} and \code{allele2} (character, numeric and numeric respectively).

Examples of these 3 data formats are objects \code{db}, \code{epg} and \code{Cgt}, respectively, in \code{test2data}.

A \bold{allele count array} is an alternative format for a genotype as a named list of vectors, one for each marker. Each vector gives the number of each allele in the genotype, with the alleles listed in the order in which they appear in the \code{data} component of the relevant mixture object. 

An \bold{IBD pattern distribution} is a list with components \code{pr} (a numerical vector) and \code{patt} (an integer matrix with \code{nrow(patt)==length(pr)} and an even number of columns). The elements of \code{pr} are the probabilities of the IBD patterns in the corresponding rows of \code{patt}. Adjacent pairs of columns encode the genotypes of different individuals; equal elements in any row determine equality of the alleles; different elements denote independent draws from the gene pool. If the component \code{pr} is missing, functions \code{rpt.IBD} and \code{rpt.typed.relatives} assume the probabilities are equal.
}

\examples{
require(KinMixLite)
require(ribd)
data(test2data)
data(NGMDyes)

C<-50

## Fit 2-person mixture - baseline model

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.7,U2=0.3)))
baseline<-logL(mixD)(pars)

## Fit 2-person mixture model in which contributor 1 is parent of a typed individual Cgt

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,'parent',list(c=Cgt),targets=c('f','c'),contrib='f') 
log10LR<-(logL(mixD)(pars)-baseline)/log(10)
cat('log10 LR',log10LR,'\n')

## Fit 3-person mixture, where 1st and 3rd contributors are siblings

mixD<-DNAmixture(list(epg),k=3,C=list(C),database=db)
pars<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.6,U2=0.3,U3=0.1)))
baseline<-logL(mixD)(pars)

mixD<-DNAmixture(list(epg),k=3,C=list(C),database=db,triangulate=FALSE,compile=FALSE)
rpt.IBD(mixD,'sibs',targets=c('b1','b2'),contribs=c('b1','u','b2')) 
log10LR<-(protected(logL(mixD)(pars))-baseline)/log(10)
cat('log10 LR',log10LR,'\n')
}
