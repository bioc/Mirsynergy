\name{tcga_brca_testdata}
\alias{tcga_brca_testdata}
\docType{data}
\title{
Breast cancer expression test data from TCGA
}
\description{
Test data of 2661 mRNA and 142 miRNA (i.e., 0.2 of the whole data) across 15
tumor samples from breast cancer (BRCA) patients. The data were downloaded from
TCGA (The Cancer Genome Atlas). The full BRCA expression data contain expression
measurements for 13306 and 710 distinct mRNAs and miRNAs across 331 samples.
}
\format{
  A list containing the follow items:
  \describe{
    \item{X,Z}{N-by-T and M-by-T Expression matrices for N mRNA and M miRNA
measured across T samples}
    \item{C,H}{N-by-M and N-by-N matrices for sequence-based miRNA-targets
downloaded from TargetScanHuman 6.2 and gene-gene interactions from TRANSFAC and
BioGrid invoving transcription factor binding sites (TFBS) and protein-protein
interactions, respectively.}
  }
}
\references{
Cancer Genome Atlas Research Network (2008). Comprehensive genomic
characterization defines human glioblastoma genes and core pathways. Nature,
455(7216), 1061-1068.

Friedman, R. C., Farh, K. K.-H., Burge, C. B., and Bartel, D. P. (2009). Most
mammalian mRNAs are conserved targets of microRNAs. Genome Research, 19(1),
92-105.
}
\examples{
load(system.file("extdata/tcga_brca_testdata.RData", package="Mirsynergy"))
}
\keyword{datasets}
