\name{mirsynergy}
\alias{mirsynergy}
\title{
Detect synergistic miRNA regulatory modules by overlapping neighbourhood
expansion
}
\description{
Detect synergistic miRNA regulatory modules by overlapping neighbourhood
expansion using a deterministic overlapping clustering algorithm adapted from a
recently developed framework. Mirsynergy operates in two stages that first forms
MRM based on co-occurring miRNAs and then expand the MRM by greedily including
(excluding) mRNA into (from) the MRM to maximize the synergy score, which is a
function of miRNA-mRNA and gene-gene interactions.
}
\usage{
mirsynergy(W, H, alpha = 2, merge.tol = 0.8, 
  density1.tol = 1e-2, density2.tol=5e-3, verbose = FALSE)
}
\arguments{
  \item{W}{
An N by M edge weight matrix containing interaction strength between N mRNA and
M miRNA.
}
  \item{H}{
An N by N edge weight matrix containing the binary interaction among the N mRNA
(genes).
}
  \item{alpha}{
Penalty for including a node into the growing module (advanced option). See
manuscript or Nepusz et al. (2012) for more details.
}
  \item{merge.tol}{
Threshold with range [0,1] to merge modules based on the percentage of nodes
shared between the two modules.
}
  \item{density1.tol}{
Threshold with range [0,1] to filter modules by the density function at stage 1
clustering.
}
  \item{density2.tol}{
Threshold with range [0,1] to filter modules by the density function 2 at stage
2 clustering.
}
  \item{verbose}{
Binary indicator to show running info.
}}

\details{
The weight matrix W can be obtained by various approaches such as Pearson
correlation or linear regression on mRNA and miNRA expression profiles across
multiple samples. Matrix H can be obtained from public database such as TRANSFAC
and BioGrid.
}
\value{
A nested list containing each item as a miRNA regulatory module (MRM). Each item
itself is a list containing the following information:
\item{miRNA}{miRNA included in the MRM}
\item{mRNA}{mRNA included in the MRM}
\item{v.in}{miRNA and mRNA}
\item{v.bound}{miRNA and mRNA disregard or excluded from the MRM but still have
nonzero connection with the internal nodes}
\item{card.m}{Number of miRNA in the MRM}
\item{card.t}{Number of mRNA targets in the MRM}
\item{card}{Total number of miRNA and mRNA targets in the MRM}
\item{density}{Density of the MRM}
}
\references{
Nepusz, T., Yu, H., & Paccanaro, A. (2012). Detecting overlapping protein
complexes in protein-protein interaction networks. Nature Methods, 9(5),
471-472. doi:10.1038/nmeth.1938
}
\author{
Yue Li
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
##
# simulate N mRNA and M miRNA and their interaction matrices
load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))

# run mirsynergy clustering
V <- mirsynergy(W, H, verbose=TRUE)

summary_modules(V)
}
\keyword{clustering}
\keyword{microRNA}

