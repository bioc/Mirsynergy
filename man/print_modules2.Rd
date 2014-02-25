\name{print_modules2}
\alias{print_modules2}
\title{
Print basic information of the formed modules
}
\description{
Print basic information of the formed modules.
}
\usage{
print_modules2(V)
}
\arguments{
  \item{V}{
Output from \code{\link{mirsynergy}}.
}
}
\author{
Yue Li
}
\examples{
load(system.file("extdata/toy_modules.RData", package="Mirsynergy"))

print_modules2(V)

}
