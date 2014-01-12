\name{kaksAA-class}
\docType{class}
\alias{kaksAA-class}
\alias{show,kaksAA-method}
\alias{seq_num.CorMut}
\alias{kaks.CorMut}
\alias{q.CorMut}
\alias{ka.CorMut}
\alias{ks.CorMut}
% Accessor methods:
\alias{seq_num.CorMut,kaksAA-method}
\alias{kaks.CorMut,kaksAA-method}
\alias{lod.CorMut,kaksAA-method}
\alias{q.CorMut,kaksAA-method}
\alias{ka.CorMut,kaksAA-method}
\alias{ks.CorMut,kaksAA-method}
\title{Class "kaksAA"}
\description{kaksAA class for representing the kaksAA results}
\section{Slots}{
    \describe{
    \item{\code{seq_num}:}{Object of class \code{"list"}. The number of sequences}
    \item{\code{kaks}:}{Object of class \code{"list"}. kaks value}
    \item{\code{lod}:}{Object of class \code{"list"}. LOD confidence score, refer to vignette for details}
    \item{\code{q}:}{Object of class \code{"list"}. q value, refer to vignette for details}
    \item{\code{ka}:}{Object of class \code{"list"}.ka value of kaks ratio, i.e. the samples count of the nonsynonymous mutation. }
    \item{\code{ks}:}{Object of class \code{"list"}. ks value of kaks ratio, i.e. the samples count of the synonymous mutation.}
}
}

\section{Methods}{
  \describe{
    \item{filterSites}{\code{signature(object = "kaksAA")}: Filter positive selection amino mutations for objects of kaksAA class.}
  }
}
\seealso{
\code{\link{filterSites}},\code{\link{kaksAA}}
}
\examples{
#examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
#example=seqFormat(examplefile)
#result=kaksAA(example)
}
\keyword{classes}
