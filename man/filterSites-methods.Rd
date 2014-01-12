\name{filterSites-methods}
\docType{methods}
\alias{filterSites}
\alias{filterSites-methods}
\alias{filterSites.kaksCodon}
\alias{filterSites.kaksAA}
\alias{filterSites.ckaks}
\alias{filterSites.MI}
\alias{filterSites.JI}
\alias{filterSites,kaksCodon-method}
\alias{filterSites,kaksAA-method}
\alias{filterSites,ckaks-method}
\alias{filterSites,MI-method}
\alias{filterSites,JI-method}
\title{filterSites methods}
\description{
filterSites implement filterSites methods for kaksCodon,kaksAA,ckaks and MI instances respectively.The functions filter the results of corresponding objects.
}
\section{Methods}{
\describe{
\item{x = "kaksCodon"}{
\code{filterSites(x,lod_cut = 2,freq_cut=0.01)}:
 Filter the results of kaksCodon object,two parameters are provided to control the output, lod_cut can assign an cutoff for lod confidence score to filter, and freq_cut assign the cutoff frequence of mutation, thus low-frequency mutations are ignored in the later analyses.
}
\item{x = "kaksAA"}{
\code{filterSites(x,lod_cut = 2,freq_cut=0.01)}:
 Filter the results of kaksAA object,two parameters are provided to control the output, lod_cut can assign an cutoff for lod confidence score to filter, and freq_cut assign the cutoff frequence of mutation, thus low-frequency mutations are ignored in the later analyses.
}
\item{x = "ckaks"}{
\code{filterSites(x,lod_cut=2)}:
Filter the results of ckaks object,lod_cut can assign an cutoff for lod confidence score to filter.
}
\item{x = "MI"}{
\code{filterSites(x,p_cut=0.05)}:
Filter the results of MI object,p_cut can assign an p value cutoff to filter the correlated mutations.
}
\item{x = "JI"}{
\code{filterSites(x,p_cut=0.05)}:
Filter the results of MI object,p_cut can assign an p value cutoff to filter the correlated mutations.
}
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
