## load libraries
library(foreach)

#' Coefficients for half-circle FRAP, diffusion without boundary, bleached half
#'
#' @param alpha Zeros for evaluation
#' @param n Order parameter for Bessel functions
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @return Coefficients
bleached_no_boundary_half_coeff <- function(alpha, n, rc, rl) {
  res <- 0
  if(n==0) {res <- 2/rl^2 * besselJ(alpha*rc,1)^2/alpha^2/besselJ(alpha*rl,0)^2}
  else if(n%%2!=0) {
    m <- abs(n)
    res <- 8/pi^2 * rc^2*alpha^2*(alpha*rc/2)^(2*m) * hypergeo(m=m,z=-alpha^2*rc^2/4)^2
    res <- res / (alpha^2-m^2/rl^2)/m^4/(2+m)^2/gamma(m)^2/besselJ(alpha*rl,m)^2/rl^2
  }
  return(res)
}

#' Coefficients for half-circle FRAP, diffusion without boundary, non-bleached half
#'
#' @param alpha Zeros for evaluation
#' @param n Order parameter for Bessel functions
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @return Coefficients
nonbleached_no_boundary_half_coeff <- function(alpha, n, rc, rl) {
  res <- 0
  if(n==0) {res <- 2/rl^2 * besselJ(alpha*rc,1)^2/alpha^2/besselJ(alpha*rl,0)^2}
  else if(n%%2!=0) {
    m <- abs(n)
    res <- -8/pi^2 * rc^2*alpha^2*(alpha*rc/2)^(2*m) * hypergeo(m=m,z=-alpha^2*rc^2/4)^2
    res <- res / (alpha^2-m^2/rl^2)/m^4/(2+m)^2/gamma(m)^2/besselJ(alpha*rl,m)^2/rl^2
  }
  return(res)
}

#' Half-circle FRAP function, diffusion without boundary, bleached half
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @return FRAP curve
#' @export
bleached_no_boundary_half_d_fit <- function(t) {
  function(d, rc, rl) {
    amax <- 50
    nmax <- 30
    prec <- 1e-3

    cnt <- 1
    res <- matrix(0, nrow=2*nmax+1, ncol=length(t))
    for(n in -nmax:nmax) {
      # find zeros for a given n
      alpha <- zeros(rl, 0, n, amax, prec)
      alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

      # calculate function
      s <- bleached_no_boundary_half_coeff(alpha, n, rc, rl)
      res[cnt,] <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
        sum(s*exp(-alpha^2*d*t[i]))
      }

      # increase counter
      cnt <- cnt+1
    }
    return(1 - rc^2/2/rl^2 - colSums(res, na.rm = T))
  }
}

#' Half-circle FRAP function, diffusion without boundary, non-bleached half
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @return FRAP curve
#' @export
nonbleached_no_boundary_half_d_fit <- function(t) {
  function(d, rc, rl) {
    amax <- 50
    nmax <- 30
    prec <- 1e-3

    cnt <- 1
    res <- matrix(0, nrow=2*nmax+1, ncol=length(t))
    for(n in -nmax:nmax) {
      # find zeros
      alpha <- zeros(rl, 0, n, amax, prec)
      alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

      # calculate function
      s <- nonbleached_no_boundary_half_coeff(alpha, n, rc, rl)
      res[cnt,] <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
        sum(s*exp(-alpha^2*d*t[i]))
      }

      # increase counter
      cnt <- cnt+1
    }
    return(1 - rc^2/2/rl^2 - colSums(res, na.rm = T))
  }
}

#' Half-circle FRAP function, reaction-diffusion without boundary, bleached half
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @param kon Pseudo-association rate
#' @param koff Dissociation rate
#' @return FRAP curve
#' @export
bleached_no_boundary_half_rd_fit <- function(t) {
  function(d, rc, rl, kon, koff) {
    amax <- 50
    nmax <- 30
    prec <- 1e-3

    cnt <- 1
    res <- matrix(0, nrow=2*nmax+1, ncol=length(t))
    for(n in -nmax:nmax) {
      # find zeros for a given n
      alpha <- zeros(rl, 0, n, amax, prec)
      alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

      # define quantities for fit function
      feq <- koff/(kon+koff)
      v <- sqrt(1/4*(d*alpha^2+kon+koff)^2-koff*d*alpha^2)
      w <- 1/2*(d*alpha^2+kon+koff)
      a <- (koff+kon-v-w)*(v-w)/2/koff/v
      b <- (koff+kon+v-w)*(v+w)/2/koff/v

      # calculate function
      s <- bleached_no_boundary_half_coeff(alpha, n, rc, rl)
      res[cnt,] <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
        sum(s*(a*exp(-(w+v)*t[i])+b*exp(-(w-v)*t[i])))
      }

      # increase counter
      cnt <- cnt+1
    }
    return(1 - rc^2/2/rl^2 - feq*colSums(res, na.rm = T))
  }
}

#' Half-circle FRAP function, reaction-diffusion without boundary, non-bleached half
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @param kon Pseudo-association rate
#' @param koff Dissociation rate
#' @return FRAP curve
#' @export
nonbleached_no_boundary_half_rd_fit <- function(t) {
  function(d, rc, rl, kon, koff) {
    amax <- 50
    nmax <- 30
    prec <- 1e-3

    cnt <- 1
    res <- matrix(0, nrow=2*nmax+1, ncol=length(t))
    for(n in -nmax:nmax) {
      # find zeros
      alpha <- zeros(rl, 0, n, amax, prec)
      alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

      # define quantities for fit function
      feq <- koff/(kon+koff)
      v <- sqrt(1/4*(d*alpha^2+kon+koff)^2-koff*d*alpha^2)
      w <- 1/2*(d*alpha^2+kon+koff)
      a <- (koff+kon-v-w)*(v-w)/2/koff/v
      b <- (koff+kon+v-w)*(v+w)/2/koff/v

      # calculate function
      s <- nonbleached_no_boundary_half_coeff(alpha, n, rc, rl)
      res[cnt,] <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
        sum(s*(a*exp(-(w+v)*t[i])+b*exp(-(w-v)*t[i])))
      }

      # increase counter
      cnt <- cnt+1
    }
    return(1 - rc^2/2/rl^2 - feq*colSums(res, na.rm = T))
  }
}
