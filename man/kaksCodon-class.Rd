\name{kaksCodon-class}
\docType{class}
\alias{kaksCodon-class}
\alias{show,kaksCodon-method}
% Accessor methods:
\alias{seq_num.CorMut,kaksCodon-method}
\alias{kaks.CorMut,kaksCodon-method}
\alias{lod.CorMut,kaksCodon-method}
\alias{q.CorMut,kaksCodon-method}
\alias{ka.CorMut,kaksCodon-method}
\alias{ks.CorMut,kaksCodon-method}
\title{Class "kaksCodon"}
\description{kaksCodon class for representing the kaksCodon results}
\section{Slots}{
    \describe{
    \item{\code{seq_num}:}{Object of class \code{"list"}. The number of sequences}
    \item{\code{kaks}:}{Object of class \code{"numeric"}. kaks value}
    \item{\code{lod}:}{Object of class \code{"numeric"}. LOD confidence score, refer to vignette for details}
    \item{\code{q}:}{Object of class \code{"numeric"}. q value, refer to vignette for details}
    \item{\code{ka}:}{Object of class \code{"numeric"}.ka value of kaks ratio, i.e. the samples count of the nonsynonymous mutation. }
    \item{\code{ks}:}{Object of class \code{"numeric"}. ks value of kaks ratio, i.e. the samples count of the synonymous mutation.}
}
}

\section{Methods}{
  \describe{
    \item{filterSites}{\code{signature(object = "kaksCodon")}: Filter positive selection sites for objects of kaksCodon class.}
  }
}
\seealso{\code{\link{kaksCodon}},\code{\link{filterSites}}
}
\examples{
#examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
#example=seqFormat(examplefile)
#result=kaksCodon(example)
}
\keyword{classes}
