## ----echo=FALSE---------------------------------------------------------------
library(knitr)
library(qsimulatR)
knitr::opts_chunk$set(fig.align='center',
                      comment='')

## -----------------------------------------------------------------------------
N <- 3
v <- seq(1:2^N)
v <- v/sqrt(sum(v^2))

w <- fft(v, inverse=TRUE)/sqrt(length(v))

## -----------------------------------------------------------------------------
x <- qstate(N, coefs=as.complex(v))
x <- H(3) * x
x <- cqgate(bits=c(2, 3), gate=S(3)) * x
x <- cqgate(bits=c(1, 3), gate=Tgate(3)) * x
x <- H(2) * x
x <- cqgate(bits=c(1, 2), gate=S(2)) * x
x <- H(1) * x
x <- SWAP(c(1,3)) * x

## -----------------------------------------------------------------------------
plot(x)

## -----------------------------------------------------------------------------
sum(x@coefs*Conj(x@coefs))

## -----------------------------------------------------------------------------
sqrt(sum((w - x@coefs)*Conj(w - x@coefs)))

## -----------------------------------------------------------------------------
Tdagger <- function(bit) {
  return(methods::new("sqgate",
                      bit=as.integer(bit),
                      M=array(as.complex(c(1., 0, 0, exp(-1i*pi/4))), dim=c(2,2)),
                      type="Tdag"))
}
Sdagger <- function(bit) {
  return(methods::new("sqgate",
                      bit=as.integer(bit),
                      M=array(as.complex(c(1,0,0,-1i)), dim=c(2,2)),
                      type="Sdag"))
}

## -----------------------------------------------------------------------------
z <- qstate(N, coefs=x@coefs)
z <- SWAP(c(1,3)) * z
z <- H(1) * z
z <- cqgate(bits=c(1, 2), gate=Sdagger(2)) * z
z <- H(2) * z
z <- cqgate(bits=c(1, 3), gate=Tdagger(2)) * z
z <- cqgate(bits=c(2, 3), gate=Sdagger(2)) * z
z <- H(3) * z

## -----------------------------------------------------------------------------
plot(z)

## -----------------------------------------------------------------------------
sqrt(sum((v - z@coefs)*Conj(v - z@coefs)))

## -----------------------------------------------------------------------------
y <- qstate(N, coefs=x@coefs)
y <- H(3) * y
y <- cqgate(bits=c(2, 3), gate=Sdagger(3)) * y
y <- cqgate(bits=c(1, 3), gate=Tdagger(3)) * y
y <- H(2) * y
y <- cqgate(bits=c(1, 2), gate=Sdagger(2)) * y
y <- H(1) * y
y <- SWAP(c(1,3)) * y

plot(y)

## -----------------------------------------------------------------------------
sqrt(sum((v - y@coefs)*Conj(v - y@coefs)))

## -----------------------------------------------------------------------------
Ri <- function(bit, i, sign=+1) {
  type <- paste0("R", i)
  if(sign < 0) {
    type <- paste0("R", i, "dag")
  }
  return(methods::new("sqgate",
                      bit=as.integer(bit),
                      M=array(as.complex(c(1,0,0,exp(sign*2*pi*1i/2^i))),
                              dim=c(2,2)), type=type))
}

## ---- eval=FALSE--------------------------------------------------------------
#  qft <- function(x, inverse=FALSE) {
#    n <- x@nbits
#    y <- x
#    sign <- +1
#    if(inverse) sign <- -1
#    for(bit in c(n:1)) {
#      y <- H(bit) * y
#      if(bit > 1) {
#        for(i in c((bit-1):1)) {
#          y <- cqgate(bits=c(i, bit), gate=Ri(bit, bit-(i-1), sign=sign)) * y
#        }
#      }
#    }
#    ## reverse order
#    for(k in c(1:floor(n/2))) {
#      y <- SWAP(c(k, n-(k-1))) * y
#    }
#    return(invisible(y))
#  }

## -----------------------------------------------------------------------------
y <- qstate(N, coefs=x@coefs)
y <- qft(y, inverse=TRUE)
plot(y)
sqrt(sum((v - y@coefs)*Conj(v - y@coefs)))

