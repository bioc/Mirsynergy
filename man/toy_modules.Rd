\name{toy_modules}
\alias{toy_modules}
\docType{data}
\title{
Test data of 20 mRNA and 20 miRNA.
}
\description{
Test data of 20 mRNA and 20 miRNA generated from simulation.
}
\format{
  A list containing the follow items:
  \describe{
    \item{V}{Outputs from \code{\link{mirsynergy}}}
    \item{W,H}{Inputs used to generate V}
  }
}

\examples{
load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))
}
\keyword{datasets}
