\name{make.findings}
\alias{make.findings}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coding additional findings}
\description{Returns a function, by convention called extra.findings, that will be called 
in logL.UKX whenever the likelihood of the model is calculated using the
current method, and encodes the additional findings needed to implement the 
method; the details of the model and the extra information needed are held in the list z
}
\usage{
make.findings(z)
}
\arguments{
  \item{z}{A list specifying the additional findings; for the format, see Details}
}
\details{
Each component of the list \code{z} is a list encoding a particular type of additional finding: the first component of this (sub-)list being a character string specifying the type of finding, and the remainder of its components being named parameters giving details of the finding. The types of finding and the valid parameters of each are as follows:
\describe{
\item{\code{Male}}{\code{ind}: index of relevant contributor: which `unknown' contributor are we modelling by amending his CPTs?}
\item{\code{Female}}{\code{ind}: index of relevant contributor}
\item{\code{Rlikd}}{\code{aca}: allele count array, \code{cgt}: character string naming genotype profile data frame, \code{evid}: character string naming list with one component for each marker, whose value is the evidence}
\item{\code{Aca}}{\code{ind}: index of relevant contributor, \code{aca}: allele count array}
\item{\code{Caca}}{\code{ind}: index of relevant contributor, \code{aca}: allele count array}
\item{\code{Denom}}{no parameters}
}
If \code{z} is NULL, then there are no additional findings.
}
\value{Function encoding the additional findings}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
cgtcaca <- structure(list(D10 = c(0, 1, 0, 1), D12 = c(0, 1, 0, 1, 0)), logGt = -3.54737989184024)
 
extra.findings<-make.findings(list(
	list('Male',ind=1),
	list('Caca',aca='cgtcaca')
	))

extra.findings
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
