\name{miCodon}
\alias{miCodon}
\title{
Compute the mutual information among codons
}
\description{
Compute the mutual information among codons, the amino mutation in a site will be treated as whole in computation.
}
\usage{
miCodon(seq_formated, kaks = TRUE, lod_cut = 2, setPosition = c(),fdr=FALSE)
}
\arguments{
  \item{seq_formated}{
Formated alignment sequence. i.e. the result after the treatment of DataFormatCorMut.
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
A object of MI class will be return. miCodon includes two slots of matrix:mi and p.value, which indicate the mutual information and p value respectively.
}
\references{
Cover, T. M., Thomas, J. A. & Wiley, J. Elements of information theory. 6, (Wiley Online Library: 1991).
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
result=miCodon(example)
}