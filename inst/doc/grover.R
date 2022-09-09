## ---- echo=FALSE--------------------------------------------------------------
library(knitr)
library(qsimulatR)
knitr::opts_chunk$set(fig.align='center',
                      comment='')

## ---- echo=FALSE--------------------------------------------------------------
x <- qstate(2)
x <- cqgate(c(1,2), sqgate(bit=2, M=array(as.complex(c(1, 0, 0, 1)), dim=c(2,2)), type="Oracle")) * x
x <- Z(2)*x
x <- cqgate(c(1,2), sqgate(bit=2, M=array(as.complex(c(1, 0, 0, 1)), dim=c(2,2)), type="Oracle")) * x
plot(x)

## ---- echo=FALSE--------------------------------------------------------------
x <- qstate(2)
x <- X(1) * (H(1) * x)
x <- CNOT(c(1,2)) * x
x <- Z(2) * x
x <- H(1) * (X(1) * (CNOT(c(1,2)) * x))
plot(x)

## -----------------------------------------------------------------------------
## oracle for n=2 and x_s=2
oracle <- function(x) {
  x <- X(1) * (CCNOT(c(1,2,3)) *(X(1) * x))
  return(x)
}

## -----------------------------------------------------------------------------
## case |00>=0
x <- oracle(qstate(3))
measure(x, 3)$value
## case |01>=1
x <- oracle(X(1)*qstate(3))
measure(x, 3)$value
## case |10>=2
x <- oracle(X(2)*qstate(3))
measure(x, 3)$value
## case |11>=3
x <- oracle(X(2)*(X(1)*qstate(3)))
measure(x, 3)$value

## -----------------------------------------------------------------------------
U <- function(x) {
  x <- oracle(x)
  x <- Z(3) * x
  x <- oracle(x)
  return(x)
}
V <- function(x) {
  for(i in c(1:2)) {
    x <- H(i) * x
  }
  x <- X(1) * (X(2) * x)
  x <- CCNOT(c(1,2,3)) * x
  x <- Z(3) * x
  x <- CCNOT(c(1,2,3)) * x
  x <- X(1) * (X(2) * x)
  for(i in c(1:2)) {
    x <- H(i) * x
  }
  return(x)
}

## ---- echo=FALSE--------------------------------------------------------------
x <- qstate(3)
x <- U(x)
x <- V(x)
plot(x)

## -----------------------------------------------------------------------------
## prepare psi
psi <- H(1) * ( H(2) * qstate(3))
## apply VU
x <- U(psi)
x <- V(x)
x

