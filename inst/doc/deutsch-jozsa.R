## ----echo=FALSE---------------------------------------------------------------
library(knitr)
library(qsimulatR)

## -----------------------------------------------------------------------------
x <- X(1) * qstate(nbits=2, basis=genComputationalBasis(2, collapse=","))
x

## -----------------------------------------------------------------------------
y <- H(2) * (H(1) * x)
y

## -----------------------------------------------------------------------------
z <- CNOT(c(2, 1)) * y
z

## -----------------------------------------------------------------------------
u <- H(2) * z
u

## -----------------------------------------------------------------------------
value <- measure(u, 2)$value

## ---- fig.align='center'------------------------------------------------------
plot(u, qubitnames=c("|y>", "|x>"), cbitnames="c")

## -----------------------------------------------------------------------------
x <- X(1) * qstate(nbits=2, basis=genComputationalBasis(2, collapse=","))
y <- H(2) * (H(1) * x)
z <- X(1) * y
z
u <- H(2) * z
u
v <- measure(u, 2)
v
plot(v$psi)
value <- v$value

## ----comment=''---------------------------------------------------------------
filename <- paste0(tempdir(), "/circuit.py")
export2qiskit(u, filename=filename)
cat(readLines(filename), sep = '\n')

## ----comment=''---------------------------------------------------------------
export2qiskit(v$psi, filename=filename, import=TRUE)
cat(readLines(filename), sep = '\n')

