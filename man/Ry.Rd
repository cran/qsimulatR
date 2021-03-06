% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sqgate.R
\name{Ry}
\alias{Ry}
\title{The Ry gate}
\usage{
Ry(bit, theta = 0)
}
\arguments{
\item{bit}{integer. The bit to which to apply the gate}

\item{theta}{numeric. angle}
}
\value{
An S4 class 'sqgate' object is returned
}
\description{
The Ry gate
}
\examples{
x <- qstate(nbits=2)
z <- Ry(1, pi/4) * x
z

}
