\name{write.db}
\alias{write.db}
\alias{read.db}
\alias{write.epg}
\alias{read.epg}
\title{Functions for translating db and epg objects from pcrsim to KinMixLite}
\description{Functions for translating database and epg objects from pcrsim to KinMixLite format, used either directly or via .csv files}
\usage{
write.epg(res,file='epg.csv')

write.db(dbname,file='db.csv',C=0)

epg<-read.epg(file='epg.csv')

db<-read.db(file='db.csv')
}
\arguments{
  \item{res}{epg as a data frame, as output from \code{simCE}}
  \item{file}{file name}
  \item{dbname}{name of a database, a string as understood by \code{getDb}}
  \item{C}{threshold, peak heights less than \code{C} are suppressed}
}
\value{The format for epg and db is as described in \code{\link{KinMixLite}}.}
\details{if file name is an empty string, the epg or db object is returned as the value of the respective function, rather than written to a file.
}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{KinMixLite}}, \code{\link{SimCE}}, \code{\link{getDb}}
}
\examples{
require(KinMixLite)
data(test2data)
data(NGMDyes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
