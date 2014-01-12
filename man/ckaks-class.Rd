\name{ckaks-class}
\docType{class}
\alias{ckaks-class}
\alias{show,ckaks-method}
\alias{ckaks}
\alias{ckaks.CorMut}
\alias{lod.CorMut}
% Accessor methods:
\alias{ckaks.CorMut,ckaks-method}
\alias{lod.CorMut,ckaks-method}
\title{Class "ckaks"}
\description{ckaks class for representing the ckaksCodon/ckaksAA results}
\section{Slots}{
    \describe{
    \item{\code{ckaks}:}{Object of class \code{"matrix"}. Conditional kaks among codons or amino mutations}
    \item{\code{lod}:}{Object of class \code{"matrix"}. LOD confidence score indicates the significance of correlated mutations}
}
}

\section{Methods}{
  \describe{
    \item{filterSites}{\code{signature(object = "ckaks")}: Filter positive selection amino mutations for objects of ckaks class.}
    \item{plot}{\code{signature(object = "ckaks")}: Plot the influence network among mutations for objects of ckaks class.}
  }
}
\seealso{\code{\link{filterSites}},\code{\link{ckaksCodon}},\code{\link{ckaksAA}},\code{\link{plot.ckaks}}
}
\examples{
#examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
#example=seqFormat(examplefile)
#result=ckaksAA(example)
}
\keyword{classes}
