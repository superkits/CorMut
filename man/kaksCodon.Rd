\name{kaksCodon}
\alias{kaksCodon}
\title{
Compute kaks for individual codon
}
\description{
Compute kaks for individual codon, the different mutations in a site will be considered as whole.
}
\usage{
kaksCodon(seq_formated)
}
\arguments{
  \item{seq_formated}{
Formated alignment sequence. i.e. the result after the treatment of DataFormatCorMut.
}
}
\details{
Ka/Ks ratio was used as an indicator of selective pressure acting on a protein-coding gene. A Ka/Ks ratio of 1 indicates neutral selection, i.e., the observed ratio of non-synonymous mutations versus synonymous mutations exactly matches the ratio expected under a random mutation model. Thus, amino acid changes are neither being selected for nor against. A Ka/Ks value of <1 indicates negative selection pressure. That is to say most amino acid changes are deleterious and are selected against, producing an imbalance in the observed mutations that favors synonymous mutations. In the condition of Ka/Ks>1 , it indicates that amino acid changes are favored, i.e., they increase the organism's fitness. This unusual condition may reflect a change in the function of a gene or a change in environmental conditions that forces the organism to adapt. For example, highly variable viruses mutations which confer resistance to new antiviral drugs might be expected to undergo positive selection in a patient population treated with these drugs, such as HIV and HCV.
}
\value{
A kaksCodon instance will be returned.
}
\references{
Chen, L., Perlina, A. & Lee, C. J. Positive selection detection in 40,000 human immunodeficiency virus (HIV) type 1 sequences automatically identifies drug resistance and positive fitness mutations in HIV protease and reverse transcriptase. Journal of virology 78, 3722-3732 (2004).
}
\author{
Zhenpeng Li
}
\note{
The reference sequence should be included in in the alignment as the first sequence. Before using kaksCodon, the alignment should be treated with  DataFormatCorMut.
}
\seealso{
\code{\link{filterSites-methods}},\code{\link{kaksAA}}
}
\examples{
examplefile=system.file("extdata","PI_treatment.aln",package="CorMut")
example=seqFormat(examplefile)
result=kaksCodon(example)
}
