\name{tabular_module}
\alias{tabular_module}
\title{
Generate tabulated module assignments and nodes for input to Cystoscape.
}
\description{
Generate tabulated module assignments and nodes for input to Cystoscape.
}
\usage{
tabular_module(V, W, H, outdir)
}
\arguments{
  \item{V}{
Output from \code{\link{mirsynergy}}.
}
  \item{W}{
An N by M edge weight matrix containing interaction strength between N mRNA and 
M miRNA.
}
  \item{H}{
An N by N edge weight matrix containing the binary interaction among the N mRNA 
(genes).
}
  \item{outdir}{
Path to save the nodes and edges files.
}

}
\value{
\item{nodes}{data.frame containing edges in each indexed module. Each edge (row)
 contain the edge weights, the edge type (MMI or GGI), and which module they 
 belong to.}
\item{edges}{data.frame containing nodes in each indexed module. Each node (row)
 contain the node name, type (mRNA/miRNA), and which module they belong to.}
}
\author{
Yue Li
}

\examples{
load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))

tabular_module(V,W,H)
}



