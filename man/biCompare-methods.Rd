\name{biCompare-methods}
\docType{methods}
\alias{biCompare}
\alias{biCompare-methods}
\alias{biCompare.biCkaksCodon}
\alias{biCompare.biCkaksAA}
\alias{biCompare.biMICodon}
\alias{biCompare.biMIAA}
\alias{biCompare,biCkaksCodon-method}
\alias{biCompare,biCkaksAA-method}
\alias{biCompare,biMICodon-method}
\alias{biCompare,biMIAA-method}
\title{
biCompare implement biCompare methods for comparison of the correlated mutations between two conditions by the means of correlated mutation network
}
\description{
Comparison of the correlated mutations between two conditions by the means of correlated mutation network, two correlated mutation network will be displayed in a plot, which represent the global mutation correlation among sites or mutations in two conditions. In the plot, blue nodes indicate the positive selection nodes in the first condition, while red nodes indicate the positive selection nodes in the second condition.
}
\section{Methods}{
\describe{
\item{x = "biCkaksCodon"}{ Compare the correlated mutations(in the unit of codon) between two conditions using condition kaks method.}
\item{x = "biCkaksAA"}{ Compare the correlated mutations(in the unit of amino) between two conditions using condition kaks method.}
\item{x = "biMICodon"}{ Compare the correlated mutations(in the unit of codon) between two conditions using mutual information method.}
\item{x = "biMIAA"}{ Compare the correlated mutations(in the unit of amino) between two conditions using mutual information method.}
}
}

\value{
A biCompare instance is returned. 
}
\author{
Zhenpeng Li
}
\seealso{
\code{\link{biCompare-class}},\code{\link{plot.biCompare}}
}
\examples{
examplefile=system.file("extdata","PI_treatment_naive.aln",package="CorMut")
examplefile02=system.file("extdata","PI_treatment.aln",package="CorMut")
seq_formated=seqFormat(examplefile)
seq_formated02=seqFormat(examplefile02)
biexample=biCkaksAA(seq_formated,seq_formated02)
result=biCompare(biexample)
}
