###########################
### IMPORTS & LIBRARIES ###
###########################

library(MASS)     # mvrnorm
library(Matrix)   # chol
library(hypergeo) # hypergeo

####################
##### SIMULATE #####
####################

simulate               <- function(simclass, ...) UseMethod("simulate", simclass)
simulate.simulateclass <- function(simclass) return(simclass)

#########################
##### SIMULATE TIME #####
#########################

simtime <- function(object, ...) UseMethod("simtime", object)
simtime.simulateclass <- function(simclass, units = "auto", digits = 0) gettime(simclass$siminfo, units, digits)

simidentity <- function(simclass, ...) {
     simclass$siminfo$starttime <- simclass$siminfo$endtime <- Sys.time()
     return(simclass)
}

#####################
##### PROCESSES #####
#####################

##### Helping function #####
add0 <- function(x) {
     if (is.matrix(x) && ncol(x) > 1 && nrow(x) > 1) cbind(0,x)
     else c(0,as.numeric(x))
}

##### dW #####
sim_dW <- function(mu, Sigma, N, s) {
     dW <- array(data = NA, dim = c(N,s,length(mu))) # (N,s,W)
     for (i in 1:N) dW[i,,]  <- mvrnorm(n = s, mu = mu, Sigma = Sigma)
     return(dW)
}

sim_cordW <- function(rho, N, s, dt = 1) {
     simdW(c(0,0), matrix(c(dt,rho*dt,rho*dt,dt),2,2),N,s) # (N,s,2)
}

sim_bm <- function(n, TT = 1, N = 1) {
     bm <- t(replicate(n = N, expr = rnorm(n = ceiling(n*TT), mean = 0, sd = 1/sqrt(n))))
     if (N == 1) bm <- cumsum(bm)
     else bm <- t(apply(bm, 1, cumsum))
     add0(bm)
}
