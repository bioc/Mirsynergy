\name{summary_modules}
\alias{summary_modules}
\title{
Return summary information of the formed modules
}
\description{
Return summary information of the formed modules, which are the outputs from 
\code{\link{mirsynergy}}
}
\usage{
summary_modules(V)
}
\arguments{
  \item{V}{
Output from \code{\link{mirsynergy}}.
}
}

\value{
\item{moduleSummaryInfo}{Summary information per module}
\item{miRNA.internal}{miRNA count distribution across modules}
\item{mRNA.internal}{mRNA count distribution across modules}
}

\author{
Yue Li
}


\examples{
load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))

summary_modules(V)
}
