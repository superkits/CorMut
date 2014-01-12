\name{JI-class}
\docType{class}
\alias{JI-class}
\alias{show,JI-method}
% Accessor methods:
\alias{ji.CorMut,JI-method}
\alias{p.value.CorMut,JI-method}
\title{Class "JI"}
\description{JI class for representing the jiAA results}
\section{Slots}{
    \describe{
    \item{\code{JI}:}{Object of class \code{"matrix"}. Mutual information among codons or amino mutations}
    \item{\code{p.value}:}{Object of class \code{"matrix"}. P value for the significance of correlated mutations}
    \item{\code{OR}:}{Object of class \code{"matrix"}. Odds Ratios for the correlated mutations}
}
}

\section{Methods}{
  \describe{
    \item{filterSites}{\code{signature(object = "MI")}: Filter positive selection amino mutations for objects of MI class.}
    \item{plot}{\code{signature(object = "MI")}: Plot the influence interaction among mutations for objects of MI class.}
  }
}
\seealso{\code{\link{miAA}},\code{\link{filterSites}}
}

\examples{
#examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
#example=seqFormat(examplefile)
#result=jiAA(example)
}
\keyword{classes}
