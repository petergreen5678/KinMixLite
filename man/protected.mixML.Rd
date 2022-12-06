\name{protected.mixML}
\alias{protected.mixML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Protect against numerical errors in maximum likelihood computation}
\description{Attempts to catch numerical errors in maximum likelihood computation, by replacing logL values by a default value instead of NaN or other failures}
\usage{
protected.mixML(mixture, pars, constraints = NULL, phi.eq = FALSE, val = NULL,
  trace = FALSE, default=-999999, order.unknowns = TRUE, ...)
}
\arguments{
  \item{default}{value of logL to be used if numerical errors are encountered}
  \item{all others}{as in mixML}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
