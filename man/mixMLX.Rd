\name{mixMLX}
\alias{mixMLX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximises mixture likelihood}
\description{Replacement for mixML in DNAmixtures that calls logLX instead of logL.
}
\usage{
mixMLX(mixture, pars, constraints = NULL, phi.eq = FALSE, 
    val = NULL, trace = FALSE, order.unknowns = TRUE, ...)}
\arguments{
  \item{mixture}{Compiled DNAmixture object.}
  \item{pars}{Parameters, in \code{mixpar} format.}
  \item{constraints}{...}
  \item{..}{All as in \code{mixML}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list containing
\describe{
\item{mle}{The (suggested) MLE.}
\item{lik}{The log of the likelihood (log e).}
}
as well as the output from the optimisation.}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{mixML}}.
}
\examples{
require(KinMix)
data(test2data)
data(NGMDyes)

# set threshold C
C<-0.001

mixD<-DNAmixture(list(epg),k=2,C=list(C),database=db,dyes=list(NGMDyes))

# add evidence individual 1 is Male 

extra.findings<-make.findings(list(list('Male',ind=1)))

# find MLE's and maximised likelihood

startpar<-mixpar(rho=list(2),eta=list(100),xi=list(0.1),phi=list(c(U1=0.9,U2=0.1)))
mlDM<-mixMLX(mixD,startpar,trace=F) 
pars<-mlDM$mle
cat('\nBaseline model maximised likelihood:',mlDM$lik,'\n')
cat('and MLEs:\n')
print(mlDM$mle)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
