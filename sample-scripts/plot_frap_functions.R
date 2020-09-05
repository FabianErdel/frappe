## load package
library(frappe)

## define plot range
time <- 10^seq(-4, 1.5, by=0.1)

## define parameters for plotted curves
dcoeff <- 0.1 # diffusion coefficient (um^2/s)
r <- 1        # radius of bleached cirlce (um)


## plot curves for the "client-in-droplet" scenario
## and for comparison the "free diffusion" scenario
#
# full-circle FRAP, diffusion with/without boundary
plot (time, circle_no_boundary_full_d_fit(time)(d=dcoeff, rc=r, rl=10*r), type="l", lty=2, ylim=c(0,1), xlab="Time (s)", ylab="Norm. Intensity", main="Confined Diffusion (droplet interface)\nfull-circle")
lines(time, circle_boundary_full_fit(time)(d=dcoeff, h=0.1, radius=r), col="blue")
lines(time, circle_boundary_full_fit(time)(d=dcoeff, h=  1, radius=r), col="green")
lines(time, circle_boundary_full_fit(time)(d=dcoeff, h= 10, radius=r), col="red")

# half-circle FRAP, bleached half, diffusion with/without boundary
plot (time, bleached_no_boundary_half_d_fit(time)(d=dcoeff, rc=r, rl=10*r), type="l", lty=2, ylim=c(0,1), xlab="Time (s)", ylab="Norm. Intensity", main="Confined Diffusion (droplet interface)\nhalf-circle, bleached half")
lines(time, bleached_boundary_half_fit(time)(d=dcoeff, h=0.1, radius=r), col="blue")
lines(time, bleached_boundary_half_fit(time)(d=dcoeff, h=  1, radius=r), col="green")
lines(time, bleached_boundary_half_fit(time)(d=dcoeff, h= 10, radius=r), col="red")

# half-circle FRAP, non-bleached half, diffusion with/without boundary
plot (time, nonbleached_no_boundary_half_d_fit(time)(d=dcoeff, rc=r, rl=10*r), type="l", lty=2, ylim=c(0.5,1), xlab="Time (s)", ylab="Norm. Intensity", main="Confined Diffusion (droplet interface)\nhalf-circle, non-bleached half")
lines(time, nonbleached_boundary_half_fit(time)(d=dcoeff, h=0.1, radius=r), col="blue")
lines(time, nonbleached_boundary_half_fit(time)(d=dcoeff, h=  1, radius=r), col="green")
lines(time, nonbleached_boundary_half_fit(time)(d=dcoeff, h= 10, radius=r), col="red")



## plot FRAP curves for the "client-in-gel" scenario
## and for comparison the "free diffusion" scenario
#
# full-circle FRAP (reaction-diffusion without boundary)
plot (time, circle_no_boundary_full_d_fit(time)(d=dcoeff, rc=r, rl=10*r), type="l", lty=2, ylim=c(0,1), xlab="Time (s)", ylab="Norm. Intensity", main="Reaction-Diffusion (binding to gel)\nfull-circle")
lines(time, circle_no_boundary_full_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=5,    koff=1), col="blue")
lines(time, circle_no_boundary_full_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=50,   koff=1), col="green")
lines(time, circle_no_boundary_full_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=500,  koff=1), col="red")
lines(time, circle_no_boundary_full_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=5000, koff=1), col="gray")

# half-circle FRAP (reaction-diffusion without boundary)
plot (time, bleached_no_boundary_half_d_fit(time)(d=dcoeff, rc=r, rl=10*r), type="l", lty=2, ylim=c(0,1), xlab="Time (s)", ylab="Norm. Intensity", main="Reaction-Diffusion (binding to gel)\nfull-circle")
lines(time, bleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=5,    koff=1), col="blue")
lines(time, bleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=50,   koff=1), col="green")
lines(time, bleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=500,  koff=1), col="red")
lines(time, bleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=5000, koff=1), col="gray")

# half-circle FRAP (reaction-diffusion without boundary)
plot (time, nonbleached_no_boundary_half_d_fit(time)(d=dcoeff, rc=r, rl=10*r), type="l", lty=2, ylim=c(0.75,1), xlab="Time (s)", ylab="Norm. Intensity", main="Reaction-Diffusion (binding to gel)\nfull-circle")
lines(time, nonbleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=5,    koff=1), col="blue")
lines(time, nonbleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=50,   koff=1), col="green")
lines(time, nonbleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=500,  koff=1), col="red")
lines(time, nonbleached_no_boundary_half_rd_fit(time)(d=dcoeff, rc=r, rl=10*r, kon=5000, koff=1), col="gray")


