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
\item{x = "ckaks"}{Plot the results of ckaks object. The paramater "layout_reset" was provided to set the layout of the network. You can refer to "layout_ function" in igraph package to get the graph layouts }
\item{x = "MI"}{Plot the results of MI object. The paramater "layout_reset" was provided to set the layout of the network. You can refer to "layout_ function" in igraph package to get the graph layouts}
\item{x = "JI"}{Plot the results of JI object. The paramater "layout_reset" was provided to set the layout of the network. You can refer to "layout_ function" in igraph package to get the graph layouts}
\item{x = "biCompare"}{ Plot the results of biCompare object. The paramater "layout_reset" was provided to set the layout of the network. You can refer to "layout_ function" in igraph package to get the graph layouts}
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
