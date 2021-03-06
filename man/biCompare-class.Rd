\name{biCompare-class}
\docType{class}
\alias{biCompare-class}
\alias{method.CorMut}
\alias{positiveSite01.CorMut}
\alias{positiveSite02.CorMut}
\alias{state_1.CorMut}
\alias{statistic_1.CorMut}
\alias{state_2.CorMut}
\alias{statistic_2.CorMut}
% Accessor methods:
\alias{method.CorMut,biCompare-method}
\alias{positiveSite01.CorMut,biCompare-method}
\alias{positiveSite02.CorMut,biCompare-method}
\alias{state_1.CorMut,biCompare-method}
\alias{statistic_1.CorMut,biCompare-method}
\alias{state_2.CorMut,biCompare-method}
\alias{statistic_2.CorMut,biCompare-method}
\title{Class "biCompare"}
\description{biCompare class for representing the biCompare results}
\section{Slots}{
    \describe{
    \item{\code{method}:}{Object of class \code{"character"}. The method to compute the correlated mutations, including ckaksCodon,ckaksAA,miCodon,miAA}
    \item{\code{state_1}:}{Object of class \code{"matrix"}. The conditional kaks(conditional selection pressure) or mutual information among codons or amino mutations.}
    \item{\code{statistic_1}:}{Object of class \code{"matrix"}. LOD confidence score or p value among codons or amino mutations for the first condition. }
    \item{\code{state_2}:}{Object of class \code{"matrix"}. The conditional kaks(conditional selection pressure) or mutual information among codons or amino mutations.}
    \item{\code{statistic_2}:}{Object of class \code{"matrix"}. LOD confidence score or p value among codons or amino mutations for the second condition. }
    \item{\code{positiveSite01}:}{Object of class \code{"ANY"}. A vector of positive selection sites for the first condition.}
    \item{\code{positiveSite02}:}{Object of class \code{"ANY"}. A vector of positive selection sites for the second condition.}
}
}

\section{Methods}{
  \describe{
    \item{plot}{\code{signature(object = "biCompare")}: plot the result of biCompare class}
  }
}
\seealso{
\code{\link{biCompare}},\code{\link{plot-methods}}
}
\examples{
#examplefile=system.file("extdata","PI_treatment_naive.aln",package="CorMut")
#examplefile02=system.file("extdata","PI_treatment.aln",package="CorMut")
#seq_formated=seqFormat(examplefile)
#seq_formated02=seqFormat(examplefile02)
#biexample=biCkaksAA(seq_formated,seq_formated02)
#result=biCompare(biexample)
}
\keyword{classes}
