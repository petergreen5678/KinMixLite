\name{rni}
\alias{rni}
\title{Random number initialiser supporting spontaneous reproduction}
\description{Random number initialiser supporting spontaneous reproduction}
\usage{
rni(seed=0)
}
\arguments{
  \item{seed}{Integer, seed}
}
\value{Non-zero seed value that can be used to reproduce run subsequently}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}
\details{
This is a convenience front end to \code{set.seed}. A non-zero value of \code{seed} is passed directly to \code{set.seed}. Given a zero value (the default), the function calls \code{Sys.time} to generate an unpredictable starting value -- but the value ultimately passed to \code{set.seed} is both output using \code{cat} and returned invisibly, so can be used for unanticipated replica runs of a simulation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rni(0)
runif(6)
rni(0)
runif(6)
rni(3456)
runif(6)
rni(3456)
runif(6)
keep<-rni(0)
print(keep)
runif(6)
rni(keep)
runif(6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
