## ----echo=FALSE---------------------------------------------------------------
library(knitr)
library(qsimulatR)
knitr::opts_chunk$set(fig.align='center',
                      comment='')

## -----------------------------------------------------------------------------
t=6

## -----------------------------------------------------------------------------
epsilon <- 1/4
## note the log in base-2
digits <- t-ceiling(log(2+1/(2*epsilon))/log(2)) 
digits

## -----------------------------------------------------------------------------
2^(-digits)

## -----------------------------------------------------------------------------
x <- S(1) * (H(1) * qstate(t+1, basis=""))

## -----------------------------------------------------------------------------
alpha <- pi*3/7
s <- sin(alpha)
c <- cos(alpha)
## note that R fills the matrix columns first
M <- array(as.complex(c(c, -s, s, c)), dim=c(2,2)) 
Uf <- sqgate(bit=1, M=M, type=paste0("Uf"))

## -----------------------------------------------------------------------------
for(i in c(2:(t+1))) {
  x <- H(i) * x
}

## -----------------------------------------------------------------------------
for(i in c(2:(t+1))) {
  x <- cqgate(bits=c(i, 1),
              gate=sqgate(bit=1,
                          M=M, type=paste0("Uf", 2^(i-2)))) * x
  M <- M %*% M
}
plot(x)

## ---- fig.width=9, fig.height=4-----------------------------------------------
x <- qft(x, inverse=TRUE, bits=c(2:(t+1)))
plot(x)

## -----------------------------------------------------------------------------
xtmp <- measure(x)
cbits <- genStateNumber(which(xtmp$value==1)-1, t+1)
phi <- sum(cbits[1:t]/2^(1:t))

cbits[1:t]
phi

## -----------------------------------------------------------------------------
phi-alpha/(2*pi)

## -----------------------------------------------------------------------------
plot(2*abs(x@coefs[seq(1,128,2)])^2, type="l",
     ylab="p", xlab="state index")

## -----------------------------------------------------------------------------
x <- (H(1) * qstate(t+1, basis=""))

## -----------------------------------------------------------------------------
for(i in c(2:(t+1))) {
  x <- H(i) * x
}
M <- array(as.complex(c(c, -s, s, c)), dim=c(2,2)) 
for(i in c(2:(t+1))) {
  x <- cqgate(bits=c(i, 1),
              gate=sqgate(bit=1,
                          M=M, type=paste0("Uf", 2^(i-2)))) * x
  M <- M %*% M
}
x <- qft(x, inverse=TRUE, bits=c(2:(t+1)))

measurephi <- function(x, t) {
  xtmp <- measure(x)
  cbits <- genStateNumber(which(xtmp$value==1)-1, t+1)
  phi <- sum(cbits[1:t]/2^(1:t))
  return(invisible(phi))
}
phi <- measurephi(x, t=t)
2*pi*phi
phi-c(+alpha, 2*pi-alpha)/2/pi

## -----------------------------------------------------------------------------
plot(abs(x@coefs)^2, type="l",
     ylab="p", xlab="state index")

## -----------------------------------------------------------------------------
phi <- c()
for(i in c(1:N)) {
  phi[i] <- measurephi(x, t)
}
hist(phi, breaks=2^t, xlim=c(0,1))
abline(v=c(alpha/2/pi, 1-alpha/2/pi), lwd=2, col="red")

