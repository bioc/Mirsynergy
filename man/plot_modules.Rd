\name{plot_modules}
\alias{plot_modules}
\title{
Plot module assignments.
}
\description{
Plot as network graphical view of the output from \code{\link{mirsynergy}} 
\code{V}. NB: small network only.
}
\usage{
plot_modules(V, W, H, legend.pos = "topright", ...)
}
\arguments{
  \item{V}{
Outputs from \code{\link{mirsynergy}}
}
  \item{W}{
An N by M edge weight matrix containing interaction strength between N mRNA and 
M miRNA.
}
  \item{H}{
An N by N edge weight matrix containing the binary interaction among the N mRNA 
(genes).
}
  \item{legend.pos}{
Specify legend position
}
  \item{\dots}{
Other parameters passed to \code{\link{plot}}
}
}
\details{
Each node will be coloured by the modules they belong to. If a node belongs to 
multiple modules, it will be coloured differently based on the exact combination 
of the module indexes. The legend displays the corresponding colours and module 
assignments.
}
\author{
Yue Li
}
\note{
Only for small network. For large network please use 
\code{\link{tabular_module}} to create a pairwise table as input to Cytoscape.
}

\examples{
load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))

plot_modules(V,W,H)

}








