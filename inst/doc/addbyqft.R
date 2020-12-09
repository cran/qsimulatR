## ----echo=FALSE---------------------------------------------------------------
library(knitr)
library(qsimulatR)
knitr::opts_chunk$set(fig.align='center',
                      comment='')

## -----------------------------------------------------------------------------
Rtheta <- function(bit, theta=0.) {
  return(methods::new("sqgate", bit=as.integer(bit),
                      M=array(as.complex(c(1, 0, 0, exp(1i*theta))),
                              dim=c(2,2)), type="Rt"))
}

## -----------------------------------------------------------------------------
addbyqft <- function(x, y) {
  n <- x@nbits
  z <- qsimulatR::qft(x)
  for(j in c(1:n)) {
    z  <- Rtheta(bit=j, theta = 2*pi*y/2^(n-j+1)) * z
  }
  z <- qft(z, inverse=TRUE)
  return(invisible(z))
}

## -----------------------------------------------------------------------------
x <- qstate(5, basis=as.character(seq(0, 2^5-1)))
x
z <- addbyqft(x, 3)
z
z <- addbyqft(z, 5)
z
z <- addbyqft(z, 30)
z

