\name{ckaksCodon}
\alias{ckaksCodon}
\title{
Compute the conditional kaks(conditional selection pressure) among codons
}
\description{
Compute the conditional kaks(conditional selection pressure) among codons, the amino mutation in a site will be treated as whole in computation.
}
\usage{
ckaksCodon(seq_formated, kaks = TRUE, lod_cut = 2, setPosition = c())
}
\arguments{
  \item{seq_formated}{
Formated alignment sequence. i.e. the result after the treatment of DataFormatCorMut.
}
  \item{kaks}{
A logical variable to indicate whether kaks is turn on or off, if kaks is TRUE, conditional kaks will be computed only among positive seelction sites, or if kaks is FALSE, conditional kaks will be computed only among all sites of sequences.
}
  \item{lod_cut}{
The LOD confidence score cutoff, the default value is 2. If lod is larger than 2, it means the positive selection of individual site or the conditional selection pressure among sites are significant.
}
  \item{setPosition}{
The positions of sequence to compute. setPosition should be a vector of interger type indicating the positions to compute.
}
}
\value{
A object of ckaks class will be return. ckaks includes two slots of matrix:ckaks and lod, which indicate the ckaks and lod confidence score respectively.
}
\references{
Chen, L., Perlina, A. & Lee, C. J. Positive selection detection in 40,000 human immunodeficiency virus (HIV) type 1 sequences automatically identifies drug resistance and positive fitness mutations in HIV protease and reverse transcriptase. Journal of virology 78,3722-3732 (2004).
}
\author{
Zhenpeng Li
}
\seealso{
\code{\link{filterSites}},\code{\link{plot.ckaks}},\code{\link{ckaksAA}}
}
\examples{
examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
example=seqFormat(examplefile)
result=ckaksCodon(example)
}
