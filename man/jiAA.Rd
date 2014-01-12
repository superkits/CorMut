\name{jiAA}
\alias{jiAA}
\title{
Compute the Jaccard index among individual amino mutations
}
\description{
Compute the Jaccard index among individual amino mutations, the amino mutations in a specific position will be considered respectively.
}
\usage{
jiAA(seq_formated, kaks = TRUE, lod_cut = 2, setPosition = c(),fdr=FALSE)
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
  \item{fdr}{
Decide whether use FDR procedure to control the p value of the computation.The default value is False.
}
}
\value{
A object of JI class will be returned. JI class includes three slots: JI, p.value and OR, which indicate the Jaccard index, p value and odds ratios respectively.
}

\author{
Zhenpeng Li
}

\seealso{
\code{\link{filterSites}},\code{\link{miAA}}
}

\examples{
examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
example=seqFormat(examplefile)
result=jiAA(example)
}
