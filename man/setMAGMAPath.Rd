\name{setMAGMAPath}
\alias{setMAGMAPath}
\title{Manually set the MAGMA executable path}
\usage{
setMAGMAPath(path)
}
\arguments{
\item{path}{Absolute path to the MAGMA executable.}
}
\value{
The normalized MAGMA path, invisibly.
}
\description{
Stores a validated absolute path to the MAGMA executable for the current R
session. This is a fallback for cases where \code{conseguiR} cannot
auto-discover MAGMA on package load.
}
\examples{
\dontrun{
setMAGMAPath("/absolute/path/to/magma")
}
}
