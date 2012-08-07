\name{biCkaksCodon-class}
\docType{class}
\alias{biCkaksCodon}
\alias{biCkaksAA}
\alias{biMICodon}
\alias{biMIAA}
\alias{biCkaksCodon-class}
\alias{biCkaksAA-class}
\alias{biMICodon-class}
\alias{biMIAA-class}
\alias{show,biCkaksCodon-method}
\alias{show,biCkaksAA-method}
\alias{show,biMICodon-method}
\alias{show,biMIAA-method}
\title{Class "biCkaksCodon", "biCkaksAA", "biMICodon" and "biMIAA"}
\description{biCkaksCodon/biCkaksAA/biMICodon/biMIAA class for representing seqFormt result of two conditions}
\section{Slots}{
    \describe{
    \item{\code{seq_formated}:}{Object of class \code{"DNAStringSet"}. The result of seqFormt for the first condition }
    \item{\code{seq_formated02}:}{Object of class \code{"DNAStringSet"}. The result of seqFormt for the second condition }
}
}

\section{Methods}{
  \describe{
    \item{biCompare}{\code{signature(object = "biCkaksCodon")}: Compare the correlated mutations(in the unit of codon) between two conditions using condition kaks method.}
    \item{biCompare}{\code{signature(object = "biCkaksAA")}: Compare the correlated mutations(in the unit of amino) between two conditions using condition kaks method.}
    \item{biCompare}{\code{signature(object = "biMICodon")}: Compare the correlated mutations(in the unit of codon) between two conditions using mutual information method.}
    \item{biCompare}{\code{signature(object = "biMIAA")}: Compare the correlated mutations(in the unit of amino) between two conditions using mutual information method.}
    }
}
\seealso{
\code{\link{biCompare}},\code{\link{plot}}
}
\examples{
#examplefile=system.file("extdata","PI_treatment_naive.aln",package="CorMut")
#examplefile02=system.file("extdata","PI_treatment.aln",package="CorMut")
#seq_formated=seqFormat(examplefile)
#seq_formated02=seqFormat(examplefile02)
#biexample=biCkaksAA(seq_formated,seq_formated02)
}
\keyword{classes}
