## load library
library(foreach)

#' Coefficients for full-circle FRAP, diffusion with boundary
#'
#' @param h Permeability
#' @param alpha Zeros for evaluation
#' @param radius Radius of the bleached circle
#' @return Coefficients
circle_boundary_coeff <- function(h, alpha, radius) {
  res <- 4*h^2 / (radius^2*alpha^2*(alpha^2+h^2))
  return(res)
}


#' Full-circle FRAP function, diffusion with boundary
#'
#' @param t Time points
#' @param d Diffusion coefficient
#' @param h Permeability
#' @param radius Radius of the bleached circle
#' @return FRAP curve
#' @export
circle_boundary_full_fit <- function(t) {
  function(d, h, radius) {
    amax <- 50
    prec <- 1e-3

    res <- rep(0, length(t))

    # find zeros for n=0
    alpha <- zeros(radius, h, 0, amax, prec)
    alpha <- alpha[alpha>prec] # exclude values close to zero to avoid rounding errors

    # calculate function
    s <- circle_boundary_coeff(h, alpha, radius)
    res <- foreach(i = seq_along(t), .combine = cbind) %dopar% {
      sum(s*exp(-alpha^2*d*t[i]))
    }

    return(1-res)
  }
}
