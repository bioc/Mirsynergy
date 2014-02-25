\name{Mirsynergy-package}
\alias{Mirsynergy-package}
\alias{Mirsynergy}
\docType{package}
\title{
Mirsynergy: detect synergistic miRNA regulatory modules by overlapping
neighbourhood expansion.
}
\description{
Mirsynergy is a deterministic overlapping clustering algorithm adapted from a
recently developed framework. Mirsynergy operates in two stages that first forms
MRM based on co-occurring miRNAs and then expand the MRM by greedily including
(excluding) mRNA into (from) the MRM to maximize the synergy score, which is a
function of miRNA-mRNA and gene-gene interactions.
}
\details{
\tabular{ll}{
Package: \tab Mirsynergy\cr
Type: \tab Package\cr
Version: \tab 0.99.2\cr
Date: \tab 2014-02-06\cr
License: \tab GPL-2\cr
}
The main function \code{\link{mirsynergy}} takes as inputs the mRNA and miRNA
interaction matrix and gene-gene interaction matrix provided by the user or
generated from existing techniques described in the manuscript (in preparation).
The function then outputs a list  with each item as a miRNA regulatory module
(MRM) containing miRNA and mRNA ID and and other results for diagnostic purpose.
}
\author{
Yue Li
Maintainer: Yue Li <yueli@cs.toronto.edu>
}
\references{
Li, Y. et al. Mirsynergy: detect synergistic miRNA regulatory modules by
overlapping neighbourhood expansion. (in preparation).

Nepusz, T., Yu, H., & Paccanaro, A. (2012). Detecting overlapping protein
complexes in protein-protein interaction networks. Nature Methods, 9(5),
471-472. doi:10.1038/nmeth.1938
}

\seealso{
\code{\link{mirsynergy}}
}
\keyword{Mirsynergy}
