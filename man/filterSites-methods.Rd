\name{filterSites-methods}
\docType{methods}
\alias{filterSites}
\alias{filterSites-methods}
\alias{filterSites.kaksCodon}
\alias{filterSites.kaksAA}
\alias{filterSites.ckaks}
\alias{filterSites.MI}
\alias{filterSites,kaksCodon-method}
\alias{filterSites,kaksAA-method}
\alias{filterSites,ckaks-method}
\alias{filterSites,MI-method}
\title{filterSites methods}
\description{
filterSites implement filterSites methods for kaksCodon,kaksAA,ckaks and MI instances respectively.The functions filter the results of corresponding objects.
}
\section{Methods}{
\describe{
\item{x = "kaksCodon"}{ Filter the results of kaksCodon object}
\item{x = "kaksAA"}{ Filter the results of kaksAA object}
\item{x = "ckaks"}{Filter the results of ckaks object }
\item{x = "MI"}{Filter the results of MI object }
}
}

\author{
Zhenpeng Li
}

\examples{
examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
example=seqFormat(examplefile)
result=kaksCodon(example)
resultfilter=filterSites(result)
}

\seealso{\code{\link{kaksCodon}}, \code{\link{kaksAA}}, \code{\link{ckaksCodon}},
\code{\link{ckaksAA}},\code{\link{miCodon}},\code{\link{miAA}},\code{\link{plot-methods}}
}
\keyword{methods}
