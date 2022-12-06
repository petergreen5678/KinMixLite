\name{protected}
\alias{protected}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Catch numerical errors and return -Inf}
\description{Attempts to catch numerical erros in evaluating the expression \code{x}, delivering a default result instead of NaN or other failures}
\usage{
protected(x,default=-Inf) 
}
\arguments{
  \item{x}{expression to be evaluated, typically the log-likelihood of a modified mixture model}
  \item{default}{value to be delivered if numerical errors are encountered}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns -Inf in case of error, otherwise the value of \code{x}}
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
