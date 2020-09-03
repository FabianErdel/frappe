## load library
library(hypergeo)

#' Hypergeometric series
#'
#' @param m Parameter
#' @param z Primary argument
#' @param cutoff Cutoff value to prevent rounding errors
#' @return Hypergeometric series
hypergeo <- function(m, z, cutoff=400) {
  res <- hypergeo::genhypergeo(U=1+m/2, L=c(1+m,2+m/2), z=z)
  res[abs(z)>cutoff] <- 0 # avoid numerical problems for abs(z)>cutoff
  return(res)
}

#' Secant function used by "zeros"-function to find roots of a*J_{n-1}(a*r) + (h*r-n)*J_{n}(a*r)/r in interval (a0, a1)
#'
#' @param r Radius
#' @param h Permeability
#' @param a0 Start of search interval
#' @param a1 End of search interval
#' @param n Order parameter for Bessel functions
#' @param prec Precision
#' @param tol Tolerance
#' @param niter Number of iterations
#' @return Secant
secant <- function(r, h, a0, a1, n, prec, tol=1e-07, niter=5000) {
  lw <- a0
  up <- a1

  for (i in 1:niter) {
    y0 <- a0*besselJ(a0*r,n-1)+(h*r-n)*besselJ(a0*r,n)/r
    y1 <- a1*besselJ(a1*r,n-1)+(h*r-n)*besselJ(a1*r,n)/r

    a2 <- a1 - y1*(a1-a0)/(y1-y0)
    if(!is.finite(a2)) {return(0)}
    if(a2<prec) {a2 <- 0}

    if(!(a1==0 & a2==0) & !(a1>up & a2>up) & !(a1<lw & a2<lw) & a2*r <= 100000) {
      y2 <- a2*besselJ(a2*r,n-1)+(h*r-n)*besselJ(a2*r,n)/r
      if(abs(y2) < tol) {return(a2)}

      a0 <- a1
      a1 <- a2
    }
    else {return(0)}
  }
  return(0)
}

#' function to find roots of J'_{n}(a*r) + h*J_{n}(a*r) = a/2*J_{n-1}(a*r) - a/2*J_{n+1}(a*r) + h*J_{n}(a*r) = J'_{n}(a*r) + h*J_{n}(a*r) = a*J_{n-1}(a*r) + (h*r-n)*J_{n}(a*r)/r in interval (0,m)
#'
#' @param r Radius
#' @param h Permeability
#' @param n Order parameter for Bessel functions
#' @param m End of interval
#' @param prec Precision
#' @return Approximated zeros
zeros <- function(r, h, n, m, prec) {
  steps <- 10*ceiling(r)
  res <- rep(0, steps*m)
  for(i in 1:(steps*m)) {
    res[i] <- secant(r=r, h=h, a0=(i-1)/steps, a1=i/steps, n=n, prec=prec)
    if(res[i]<(i-1)/steps | res[i]>i/steps) {res[i] <- 0}
  }
  return(res[res>0])
}
