\name{MI-class}
\docType{class}
\alias{MI-class}
\alias{show,MI-method}
\alias{mi}
\alias{p.value}
% Accessor methods:
\alias{mi,MI-method}
\alias{p.value,MI-method}
\title{Class "MI"}
\description{MI class for representing the miCodon/miAA results}
\section{Slots}{
    \describe{
    \item{\code{mi}:}{Object of class \code{"matrix"}. Mutual information among codons or amino mutations}
    \item{\code{p.value}:}{Object of class \code{"matrix"}. P value for the significance of correlated mutations}
}
}

\section{Methods}{
  \describe{
    \item{filterSites}{\code{signature(object = "MI")}: Filter positive selection amino mutations for objects of MI class.}
    \item{plot}{\code{signature(object = "MI")}: Plot the influence interaction among mutations for objects of MI class.}
  }
}
\seealso{\code{\link{miCodon}},\code{\link{filterSites}}
}

\examples{
#examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
#example=seqFormat(examplefile)
#result=miAA(example)
}
\keyword{classes}
