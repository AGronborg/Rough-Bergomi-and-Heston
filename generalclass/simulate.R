###########################
### IMPORTS & LIBRARIES ###
###########################

library(MASS)   # mvrnorm
library(Matrix) # For chol

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

#################
##### PATHS #####
#################

# Helping function
add0 <- function(x) {
     if (is.matrix(x) && ncol(x) > 1 && nrow(x) > 1) cbind(0,x)
     else c(0,as.numeric(x))
}

# dW
simdW <- function(mu, Sigma, N, s) {
     dW <- array(data = NA, dim = c(N,s,length(mu))) # (N,s,W)
     for (i in 1:N) dW[i,,]  <- mvrnorm(n = s, mu = mu, Sigma = Sigma)
     return(dW)
}

simcordW <- function(rho, N, s, dt = 1) {
     simdW(c(0,0), matrix(c(dt,rho*dt,rho*dt,dt),2,2),N,s) # (N,s,2)
}

# fBm
cov_fbm <- function(s, t, H) {
     0.5*(abs(t)^(2*H)+abs(s)^(2*H)-abs(t-s)^(2*H))
}

covmat_fbm <- function(H, n, TT = 1) {
     s  <- ceiling(n*TT)      # time steps total
     TT <- s/n                # maturity
     dt <- 1/n               # time step
     t  <- cumsum(rep(dt, s)) # time step sequence
     
     Sigma <- matrix(data = NA, nrow = s, ncol = s)
     
     for (i in 1:s) {
          for (j in i:s) {
               Sigma[i,j] <- Sigma[j,i] <- cov_fbm(t[i], t[j], H)
          }
     }
     return(Sigma)
}

simfbm_cholesky <- function(H, n, TT = 1, N = 1) {
     s     <- ceiling(n*TT)
     
     Sigma <- covmat_fbm(H,n,TT)
     L     <- t(chol(Sigma))
     dW   <- matrix(rnorm(s*N), nrow = N, ncol = s)
     fbm  <- t(apply(dW, 1, function(row) L%*%row))
     
     add0(fbm)
}

simfbm_eigen <- function(H, n, TT = 1, N = 1) {
     Sigma <- covmat_fbm(H,n,TT)
     fbm   <- mvrnorm(n = N, mu = rep(0,s), Sigma = Sigma)
     add0(fbm)
}

simfbm <- function(H, n, TT = 1, N = 1) simfbm_cholesky(H, n, TT, N) 

# Note
# cov_fbm(0, t, H) = 0
# cov_fbm(t, t, H) = t^(2H)
# cov_fbm(s, t, H = 0.5) = s (Brownian motion)
# cov_fbm(a*s, a*t, H) = a^(2H)*cov_fbm(s, t, H)

# Test
x <- replicate(n = 10000, simfbm(0.3,5))
apply(x, 1, mean)
apply(x, 1, var)
diag(covmat_fbm(0.3, 5))

# Test
x <- simfbm(0.3, 5, N = 10000)
apply(x, 2, mean)
apply(x, 2, var)
diag(covmat_fbm(0.3, 5))

# Test
x <- simfbm_eigen(0.3, 5, N = 10000)
apply(x, 2, mean)
apply(x, 2, var)
diag(covmat_fbm(0.3, 5))


