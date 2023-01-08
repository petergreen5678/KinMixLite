\name{checkpeaks}
\alias{checkpeaks}
\title{Check and modify database to have positive frequencies for all observed peaks/alleles}
\description{Check whether database has positive frequencies for all peaks/alleles observed in epg and genotype profiles, and optionally modify db by addition of small positive frequencies so that it does, followed by renormalisation of frequencies for each allele to sum to 1.}
\usage{
checkpeaks(x,db,fix=0.003)
}
\arguments{
  \item{x}{data frame, the epg or genotype profile; see Format section in \code{\link{KinMixLite}}.}
  \item{db}{data frame, the db; see Format section in \code{\link{KinMixLite}}.}
  \item{fix}{numeric: if positive, increment to db frequency for each identified discrepant peak}
}
\value{(possibly modified) db}
\author{
Peter Green (P.J.Green@bristol.ac.uk)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(KinMixLite)
data(test2data)

db<-checkpeaks(epg,db)
db<-checkpeaks(Cgt,db)

Xgt<-data.frame(marker=c('D10','D12'),allele1=c(8,13),allele2=c(13,10))
db<-checkpeaks(Xgt,db)
db
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
