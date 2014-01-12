\name{seqFormat-methods}
\docType{methods}
\alias{seqFormat}
\alias{seqFormat-methods}
\title{
Process the multiple sequence alignment files
}
\description{
Process the multiple sequence alignment files, the results can be used for the other functional functions. Note that the reference sequence should be included as the first sequence.
}
\usage{
seqFormat(x, format = c("clustal","fasta","mase","phylip","msf"))
}
\arguments{
  \item{x}{
multiple sequence alignment files
}
  \item{format}{
a character string specifying the format of the file: "clustal","fasta","mase","phylip","msf", refer to read.alignment function in seqinr package for details.
}
}
\value{
A vector of sequences with names is returned.
}
\note{
To ensure an exact result, the reference sequence should be included in the alignment as the first sequence.
}
\author{
Zhenpeng Li
}
\examples{
examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
example=seqFormat(examplefile)
}