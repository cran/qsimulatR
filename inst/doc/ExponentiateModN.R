## ----echo=FALSE---------------------------------------------------------------
library(knitr)
library(qsimulatR)
knitr::opts_chunk$set(fig.align='center',
                      comment='')

## -----------------------------------------------------------------------------
summands <- function(x, n, N) {
  b <- as.integer(intToBits(x))
  ret <- c()
  for(i in c(1:N)) {
    s <- 0
    for(j in c(1:N)) {
      s <- s+as.integer(b[j])*2^(i+j-2)
    }
    ret[i] <- s %% n
  }
  return(ret)
}

## -----------------------------------------------------------------------------
x <- 3
summands(3, 2^3, 3)

## -----------------------------------------------------------------------------
cRtheta <- function(bits, theta=0.) {
  cqgate(bits=bits, gate=methods::new("sqgate", bit=as.integer(bits[2]),
                                      M=array(as.complex(c(1, 0, 0, exp(1i*theta))),
                                              dim=c(2,2)), type="Rt"))
}

## -----------------------------------------------------------------------------
cadd <- function(c, bits, x, y) {
  n <- length(bits)
  z <- cqft(c=c, x=x, bits=bits)
  for(i in c(1:n)) {
    z  <- cRtheta(bits=c(c, bits[i]), theta = 2*pi*y/2^(n-i+1)) * z
  }
  z <- cqft(c=c, x=z, inverse=TRUE, bits=bits)
  return(invisible(z))
}

## -----------------------------------------------------------------------------
basis <- c()
for(i in c(0:(2^4-1))) basis[i+1] <- paste0("|", i %/% 2, ">|", i %% 2, ">")
x <- H(1)*qstate(4, basis=basis)
c <- 1
bits <- c(2:4)
z <- cadd(c=c, bits=bits, x=x, y=5)
z
z <- cadd(c=c, bits=bits, x=z, y=2)
z
z <- cadd(c=c, bits=bits, x=z, y=8)
z

## -----------------------------------------------------------------------------
mult <- function(reg1, reg2, x, y, swap=TRUE) {
  stopifnot(length(reg1) == length(reg2))
  n <- length(reg2)
  s <- summands(y, 2^n, n)
  for(i in c(1:n)) {
    x <- cadd(c=reg1[i], bits=reg2, x=x, y=s[i])
  }
  if(swap) {
    for(i in c(1:n)) {
      x <- SWAP(c(reg1[i], reg2[i])) * x
    }
  }
  return(invisible(x))
}

## -----------------------------------------------------------------------------
basis <- c()
for(i in c(0:(2^3-1))) {
  for(j in c(0:(2^3-1))) {
    basis[i*2^3+j + 1] <- paste0("|", i, ">|", j, ">")
  }
}
x <- X(2)*qstate(6, basis=basis)
x
reg1 <- c(1:3)
reg2 <- c(4:6)
z <- mult(reg1, reg2, x=x, y=3)
z <- X(5) * z
z
z <- mult(reg1, reg2, x=z, y=3)
z

## -----------------------------------------------------------------------------
eEa <- function(a, b) {
  if(a == 0) return(c(b, 0, 1))
  res <- eEa(b %% a, a)
  return(c(res[1], res[3] - (b %/% a) * res[2], res[2]))
}

moduloinverse <- function(a, n) {
  res <- eEa(a=a, b=n)
  if(res[1] != 1) stop("inverse does not exist!")
  return(res[2] %% n)
}

## -----------------------------------------------------------------------------
cis.less <- function(c, bits, x, c1, a, y) {
  ## add ancilla bit as most significant bit to bits
  b <- c(bits, a)
  n <- length(b)
  ## cadd works modulo 2^n
  z <- cadd(c=c, bits=b, x=x, y=2^n-y)
  ## 'copy' overflow bit
  z <- CNOT(c(a, c1)) * z
  ## add back, resetting ancilla a to |0>
  z <- cadd(c=c, bits=b, x=z, y=y)
  return(z)
}

## -----------------------------------------------------------------------------
basis <- c()
for(i in c(0:(2^6-1))) {
  basis[i + 1] <-
    paste0("|", i %/% 8 ,">|a=",
    (i %/% 4) %% 2, ">|c1=", (i%/%2) %% 2,
    ">|c=", i%%2, ">")
}

x <- H(1)*qstate(6, basis=basis)
z <- cadd(c=1, bits=c(4,5,6), x=x, y=5)
z
## 5 < 7 -> c1 = 1
v <- cis.less(c=1, bits=c(4,5,6), x=z, c1=2, a=3, y=7)
v
## 5 > 3 -> c1 = 0
w <- cis.less(c=1, bits=c(4,5,6), x=z, c1=2, a=3, y=3)
w
## 5 < 9 -> c1 = 1
w <- cis.less(c=1, bits=c(4,5,6), x=z, c1=2, a=3, y=9)
w

## -----------------------------------------------------------------------------
caddmodN <- function(c, bits, c1, c2, a, x, y, N) {
  stopifnot(length(a) == 1 && length(c1) == 1 &&
            length(c2) == 1 &&
            length(unique(c(c1, c2, a))) == 3)
  y <- y %% N
  ## set c1=1 if x < N
  z <- cis.less(c=c, bits=bits, x=x, c1=c1, a=a, y=N)
  ## set c2=1 if x < N - y
  z <- cis.less(c=c, bits=bits, x=z, c1=c2, a=a, y=N-y)
  
  ## if c1 and not c2, x = x + y - N
  z <- X(c2) *( CCNOT(c(c1, c2, a)) * (X(c2) * z))
  z <- cadd(c=a, bits=bits, x=z, y=y - N)
  z <- X(c2) * (CCNOT(c(c1, c2, a)) * (X(c2) * z))
  
  ## if c1 and c2 add x = x + y
  z <- CCNOT(c(c1, c2, a)) * z
  z <- cadd(c=a, bits=bits, x=z, y=y)
  z <- CCNOT(c(c1, c2, a)) * z
  
  ## reset c1,2
  z <- cis.less(c=c, bits=bits, x=z, c1=c2, a=a, y=y)
  z <- CNOT(c(c1, c2)) * z
  z <- cis.less(c=c, bits=bits, x=z, c1=c1, a=a, y=N)
  return(invisible(z))
}

## -----------------------------------------------------------------------------
basis <- c()
for(i in c(0:(2^7-1))) {
  basis[i + 1] <-
    paste0("|", i %/% 16 , ">|a=", (i %/% 8) %% 2,
           ">|c2=", (i %/% 4) %% 2,
           ">|c1=", (i%/%2) %% 2, ">|c=", i%%2, ">")
}

x <- X(1)*qstate(7, basis=basis)
x
bits <- c(5,6,7)
c <- 1
c1 <- 2
c2 <- 3
a <- 4
N <- 5
z <- caddmodN(c=c, bits=bits, c1=c1, c2=c2, a=a, x=x, y=3, N=N) # 0 + 3 mod 5
z
z <- caddmodN(c=c, bits=bits, c1=c1, c2=c2, a=a, x=z, y=1, N=N) # 3 + 1 mod 5
z
z <- caddmodN(c=c, bits=bits, c1=c1, c2=c2, a=a, x=z, y=6, N=N) # 4 + 6 mod 5
z

## -----------------------------------------------------------------------------
cmultmodN <- function(c, reg1, reg2, ancillas, x, y, N) {
  stopifnot(length(reg1) == length(reg2))
  ## need 4 ancilla registers
  stopifnot(length(ancillas) == 4 &&
            length(unique(ancillas)) == 4)
  n <- length(reg2)
  ## precompute terms in the sum
  s <- summands(y, N, n)
  ## start with |x>|0>
  for(i in c(1:n)) {
    x <- CCNOT(c(c, reg1[i], ancillas[4])) * x
    x <- caddmodN(c=ancillas[4], bits=reg2,
                  c1=ancillas[1], c2=ancillas[2],
                  a=ancillas[3],
                  x=x, y=s[i], N=N)
    x <- CCNOT(c(c, reg1[i], ancillas[4])) * x
  }
  ## now |x>|xy mod N>
  for(i in c(1:n)) {
    x <- CSWAP(c(c, reg1[i], reg2[i])) * x
  }
  ## now |xy mod N>|x>
  ## -y_inv mod N
  yinv <- N - moduloinverse(a=y, n=N)
  s <- summands(yinv, N, n)
  for(i in c(1:n)) {
    x <- CCNOT(c(c, reg1[i], ancillas[4])) * x
    x <- caddmodN(c=ancillas[4], bits=reg2,
                  c1=ancillas[1], c2=ancillas[2],
                  a=ancillas[3],
                  x=x, y=s[i], N=N)
    x <- CCNOT(c(c, reg1[i], ancillas[4])) * x
  }
  ## finally |xy mod N>|0>
  return(invisible(x))
}

## -----------------------------------------------------------------------------
basis <- c()
for(i in c(0:(2^11-1))) {
  basis[i + 1] <-
    paste0("|reg1=", i %/% (32*2^3) , ">|reg2=", (i %/% 32) %% 2^3 ,
           "|anc=", (i %/% 16) %% 2,
           (i %/% 8) %% 2, (i %/% 4) %% 2,
           (i%/%2) %% 2, ">|c=", i%%2, ">")
}
x <- CNOT(c(1,10)) * (H(1)*qstate(11, basis=basis))
x
c <- 1
ancillas <- c(2:5)
reg2 <- c(6:8)
reg1 <- c(9:11)
N <- 5
z <- cmultmodN(c=c, reg1=reg1, reg2=reg2,
               ancillas=ancillas, x=x, y=3, N=N)
z
z <- cmultmodN(c=c, reg1=reg1, reg2=reg2,
               ancillas=ancillas, x=z, y=3, N=N)
z
z <- cmultmodN(c=c, reg1=reg1, reg2=reg2,
               ancillas=ancillas, x=z, y=3, N=N)
z
z <- cmultmodN(c=c, reg1=reg1, reg2=reg2,
               ancillas=ancillas, x=z, y=3, N=N)
z

## -----------------------------------------------------------------------------
cexpomodN <- function(c, reg1, reg2, ancillas, x, y, a, N) {
  stopifnot(length(reg1) == length(reg2))
  ## need 4 ancilla registers
  stopifnot(length(ancillas) == 4 &&
            length(unique(ancillas)) == 4)
  for(i in c(1:a)) {
    x <- cmultmodN(c=c, reg1=reg1, reg2=reg2,
                   ancillas=ancillas, x=x, y=y, N=N)
  }
  return(invisible(x))
}

## -----------------------------------------------------------------------------
x <- CNOT(c(1,10)) * (H(1)*qstate(11, basis=basis))
x
x <- cexpomodN(c, reg1, reg2, ancillas, x, y=3, a=4, N)
x

## -----------------------------------------------------------------------------
cexpomodN2 <- function(c, reg1, reg2, ancillas, x, y, a, N) {
  stopifnot(length(reg1) == length(reg2))
  ## need 4 ancilla registers
  stopifnot(length(ancillas) == 4 &&
            length(unique(ancillas)) == 4)
  ab <- as.integer(intToBits(a))
  n <- max(which(ab == 1))
  y2 <- y %% N
  for(i in c(1:n)) {
    if(ab[i] == 1) {
      x <- cmultmodN(c=c, reg1=reg1, reg2=reg2,
                     ancillas=ancillas, x=x, y=y2, N=N)
    }
    y2 <- ((y2%%N) * (y2%%N)) %% N # y2=y^(2^i) mod N
  }
  return(invisible(x))
}

## -----------------------------------------------------------------------------
x <- CNOT(c(1,10)) * (H(1)*qstate(11, basis=basis))
x
x <- cexpomodN2(c, reg1, reg2, ancillas, x, y=3, a=4, N)
x

