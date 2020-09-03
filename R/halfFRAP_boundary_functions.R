## load libraries
library(foreach)

#' Coefficients for half-circle FRAP, diffusion with boundary, bleached half
#'
#' @param h Permeability
#' @param alpha Zeros for evaluation
#' @param n Order parameter for Bessel functions
#' @param radius Radius of the bleached circle
#' @return Coefficients
bleached_boundary_coeff <- function(h, alpha, n, radius) {
  res <- 0
  if(n==0) {res <- 2*h^2 / (radius^2*alpha^2*(alpha^2+h^2))}
  else if(n%%2!=0) {
    m <- abs(n)
    res <- 8/pi^2 * alpha^2*(alpha*radius/2)^(2*m) * hypergeo(m=m,z=-alpha^2*radius^2/4)^2
    res <- res / (alpha^2+h^2-m^2/radius^2)/m^4/(2+m)^2/gamma(m)^2/besselJ(alpha*radius,m)^2
  }
  return(res)
}

#' Coefficients for half-circle FRAP, diffusion with boundary, non-bleached half
#'
#' @param h Permeability
#' @param alpha Zeros for evaluation
#' @param n Order parameter for Bessel functions
#' @param radius Radius of the bleached circle
#' @return Coefficients
nonbleached_boundary_coeff <- function(h, alpha, n, radius) {
  res <- 0
  if(n==0) {res <- 2*h^2 / (radius^2*alpha^2*(alpha^2+h^2))}
  else if(n%%2!=0) {
    m <- abs(n)
    res <- -8/pi^2 * alpha^2*(alpha*radius/2)^(2*m) * hypergeo(m=m,z=-alpha^2*radius^2/4)^2
    res <- res / (alpha^2+h^2-m^2/radius^2)/m^4/(2+m)^2/gamma(m)^2/besselJ(alpha*radius,m)^2
  }
  return(res)
}

#' Half-circle FRAP function, diffusion with boundary, bleached half
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param h Permeability
#' @param radius Radius of the bleached circle
#' @return FRAP curve
#' @export
bleached_boundary_half_fit <- function(t) {
  function(d, h, radius) {
    amax <- 50
    nmax <- 30
    prec <- 1e-3

    cnt <- 1
    res <- matrix(0, nrow=2*nmax+1, ncol=length(t))
    for(n in -nmax:nmax) {
      # find zeros for a given n
      alpha <- zeros(radius, h, n, amax, prec)
      alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

      # calculate function
      s <- bleached_boundary_coeff(h, alpha, n, radius)
      res[cnt,] <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
        sum(s*exp(-alpha^2*d*t[i]))
      }

      # increase counter
      cnt <- cnt+1
    }
    return(1-colSums(res, na.rm = T))
  }
}


#' Half-circle FRAP function, diffusion with boundary, non-bleached half
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param h Permeability
#' @param radius Radius of the bleached circle
#' @return FRAP curve
#' @export
nonbleached_boundary_half_fit <- function(t) {
  function(d, h, radius) {
    amax <- 50
    nmax <- 30
    prec <- 1e-3

    cnt <- 1
    res <- matrix(0, nrow=2*nmax+1, ncol=length(t))
    for(n in -nmax:nmax) {
      # find zeros
      alpha <- zeros(radius, h, n, amax, prec)
      alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

      # calculate function
      s <- nonbleached_boundary_coeff(h, alpha, n, radius)
      res[cnt,] <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
        sum(s*exp(-alpha^2*d*t[i]))
      }

      # increase counter
      cnt <- cnt+1
    }
    return(1-colSums(res, na.rm = T))
  }
}
