## load library
library(foreach)

#' Full-circle FRAP function, diffusion
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param radius Radius of the bleached circle
#' @return FRAP curve
#' @export
soumpasis <- function(t) {
  function(d, radius) {
    tau <- radius^2/d
    res <- exp(-tau/2/t)*(besselI(tau/2/t,0)+besselI(tau/2/t,1))
    res[t==0] <- 0 # avoid rounding errors for t=0
    return(res)
  }
}

#' Full-circle FRAP function, reaction-dominant
#'
#' @param t Time points
#' @param kon Pseudo-association rate
#' @param koff Dissociation rate
#' @return FRAP curve
#' @export
reaction <- function(t) {
  function(kon, koff) {
    res <- 1-kon/(kon+koff)*exp(-koff*t)
    return(res)
  }
}

#' Coefficients for full-circle FRAP, diffusion without boundary
#'
#' @param alpha Zeros for evaluation
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @return Coefficients
circle_no_boundary_full_coeff <- function(alpha, rc, rl) {
  res <- 4*besselJ(alpha*rc,1)^2/(alpha^2*rl^2*besselJ(alpha*rl,0)^2)
  return(res)
}

#' Full-circle FRAP function, diffusion without boundary
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @return FRAP curve
#' @export
circle_no_boundary_full_d_fit <- function(t) {
  function(d, rc, rl) {
    amax <- 50
    prec <- 1e-3

    res <- rep(0, length(t))

    # find positive zeros for n=0
    alpha <- zeros(rl, 0, 0, amax, prec)
    alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

    # calculate function
    s <- circle_no_boundary_full_coeff(alpha, rc, rl)
    res <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
      sum(s*exp(-alpha^2*d*t[i]))
    }

    return(1-rc^2/rl^2-res)
  }
}

#' Full-circle FRAP function, reaction-diffusion without boundary
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param rc Radius of the bleached circle
#' @param rl Radius of the entire reservoir
#' @param kon Pseudo-association rate
#' @param koff Dissociation rate
#' @return FRAP curve
#' @export
circle_no_boundary_full_rd_fit <- function(t) {
  function(d, rc, rl, kon, koff) {
    amax <- 50
    prec <- 1e-3

    res <- rep(0, length(t))

    # find positive zeros for n=0
    alpha <- zeros(rl, 0, 0, amax, prec)
    alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

    # define quantities for fit function
    feq <- koff/(kon+koff)
    v <- sqrt(1/4*(d*alpha^2+kon+koff)^2-koff*d*alpha^2)
    w <- 1/2*(d*alpha^2+kon+koff)
    a <- (koff+kon-v-w)*(v-w)/2/koff/v
    b <- (koff+kon+v-w)*(v+w)/2/koff/v

    # calculate function
    s <- circle_no_boundary_full_coeff(alpha, rc, rl)
    res <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
      sum(s*(a*exp(-(w+v)*t[i])+b*exp(-(w-v)*t[i])))
    }

    return(1-rc^2/rl^2-feq*res)
  }
}

