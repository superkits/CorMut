\name{plot}

\alias{plot}
\alias{plot.ckaks}
\alias{plot.MI}
\alias{plot.JI}
\alias{plot.biCompare}
\alias{plot,ckaks-method}
\alias{plot,MI-method}
\alias{plot,JI-method}
\alias{plot,biCompare-method}
\title{ Plot methods for CorMut package}
\description{
Plot implement plot methods for ckaksCodon,ckaksAA,miCodon, miAA and biCompare objects respectively. Plot visualize the mutation correlation among sites or amino mutations.
}
\section{Methods}{
\describe{
\item{x = "ckaks"}{Plot the results of ckaks object }
\item{x = "MI"}{Plot the results of MI object }
\item{x = "JI"}{Plot the results of JI object }
\item{x = "biCompare"}{ Plot the results of biCompare object}
}
}

\author{
Zhenpeng Li
}

\examples{
examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
example=seqFormat(examplefile)
result=ckaksAA(example)
plot(result)
}
\seealso{
\code{\link{ckaksCodon}},\code{\link{ckaksAA}},\code{\link{miCodon}},\code{\link{miAA}},\code{\link{jiAA}}
}
\keyword{methods}
