\name{miAA}
\alias{miAA}
\title{
Compute the mutual information among individual amino mutations
}
\description{
Compute the mutual information among individual amino mutations, the amino mutations in a specific position will be considered respectively.
}
\usage{
miAA(seq_formated, kaks = TRUE, lod_cut = 2, setPosition = c())
}
\arguments{
  \item{seq_formated}{
Formated multiple alignment sequence. i.e. the result after the treatment of DataFormatCorMut.
}
  \item{kaks}{
A logical variable to indicate whether kaks is turn on or off, if kaks is TRUE, conditional kaks will be computed only among positive seelction sites, or if kaks is FALSE, conditional kaks will be computed only among all sites of sequence.
}
  \item{lod_cut}{
The LOD confidence score cutoff, the default value is 2. If lod is larger than 2, it means the positive selection of individual site or the conditional selection pressure among sites are significant.
}
  \item{setPosition}{
The positions of sequence to compute. setPosition should be a vector of interger type.
}
}
\value{
A object of MI class will be returned. miAA includes two slots of matrix:mi and p.value, which indicate the mutual information and p value respectively.
}
\references{
Cover, T. M., Thomas, J. A. & Wiley, J. Elements of information theory. 6, (Wiley Online Library: 1991).
}
\author{
Zhenpeng Li
}

\seealso{
\code{\link{filterSites}},\code{\link{miCodon}}
}

\examples{
examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
example=seqFormat(examplefile)
result=miAA(example)
}
