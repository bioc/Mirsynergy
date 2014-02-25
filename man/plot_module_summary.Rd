\name{plot_module_summary}
\alias{plot_module_summary}
\title{
Plot module statistics
}
\description{
Plot module statistics using output from \code{\link{mirsynergy}}.
}
\usage{
plot_module_summary(V)
}
\arguments{
  \item{V}{
Outputs from \code{\link{mirsynergy}}
}
}
\author{
Yue Li
}

\examples{
load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))

plot_module_summary(V)
}
