###########################
### IMPORTS & LIBRARIES ###
###########################

library(MASS)   # mvrnorm
library(Matrix) # chol
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

#####################
##### PROCESSES #####
#####################

##### Helping function #####
add0 <- function(x) {
     if (is.matrix(x) && ncol(x) > 1 && nrow(x) > 1) cbind(0,x)
     else c(0,as.numeric(x))
}

##### dW #####
simdW <- function(mu, Sigma, N, s) {
     dW <- array(data = NA, dim = c(N,s,length(mu))) # (N,s,W)
     for (i in 1:N) dW[i,,]  <- mvrnorm(n = s, mu = mu, Sigma = Sigma)
     return(dW)
}

simcordW <- function(rho, N, s, dt = 1) {
     simdW(c(0,0), matrix(c(dt,rho*dt,rho*dt,dt),2,2),N,s) # (N,s,2)
}

##### fBm #####
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
     s     <- ceiling(n*TT)
     Sigma <- covmat_fbm(H,n,TT)
     fbm   <- mvrnorm(n = N, mu = rep(0,s), Sigma = Sigma)
     add0(fbm)
}

simfbm <- function(H, n, TT = 1, N = 1) simfbm_cholesky(H, n, TT, N) 

##### VOLTERRA #####
G <- function(x, H) {
     gamma <- 0.5 - H
     2*H*integrate(f = function(s) 1/((1-s)^gamma*(x-s)^gamma), lower = 0, upper = 1)$value
     # slower alternative:
     # Re( (1-2*gamma)/(1-gamma) * x^(-gamma) * hypergeo(1, gamma, 2 - gamma, 1/x) )
}

cov_volterra <- function(s, t, H) {
     min(s,t)^(2*H)*G(max(t/s,s/t), H)
}
cov_bm <- function(s, t) min(s,t)

cov_bm_volterra <- function(bm_t, volterra_t, H, rho) {
     DH <- sqrt(2*H)/(H+0.5)
     rho*DH*(volterra_t^(H+0.5)-(volterra_t-min(volterra_t-bm_t))^(H+0.5))
}

covmat_volterra <- function(H, n, TT = 1) {
     s  <- ceiling(n*TT)      # time steps total
     TT <- s/n                # maturity
     dt <- 1/n                # time step
     t  <- cumsum(rep(dt, s)) # time step sequence
     
     Sigma <- matrix(data = NA, nrow = s, ncol = s)
     
     for (i in 1:s) {
          for (j in i:s) {
               Sigma[i,j] <- Sigma[j,i] <- cov_volterra(t[i], t[j], H)
          }
     }

     return(Sigma)
}

covmat_bm_volterra <- function(H, rho, n, TT = 1) {
     s  <- ceiling(n*TT)      # time steps total
     TT <- s/n                # maturity
     dt <- 1/n                # time step
     t  <- cumsum(rep(dt, s)) # time step sequence
     
     Sigma <- matrix(data = NA, nrow = 2*s, ncol = 2*s)
     
     for (i in 1:s) {
          for (j in i:s) {
               Sigma[i,j] <- Sigma[j,i] <- cov_bm(t[i], t[j])
          }
     }
     for (i in 1:s) {
          for (j in 1:s) {
               Sigma[i,s+j] <- Sigma[s+j,i] <- cov_bm_volterra(t[i], t[j], H, rho)
          }
     }
     for (i in 1:s) {
          for (j in i:s) {
               Sigma[s+i,s+j] <- Sigma[s+j,s+i] <- cov_volterra(t[i], t[j], H)
          }
     }
     return(Sigma)
}

sim_volterra_cholesky    <- function(H, n, TT = 1, N = 1) {
     s     <- ceiling(n*TT)
     
     Sigma <- covmat_volterra(H, n, TT)
     L     <- t(chol(Sigma))
     dW    <- matrix(rnorm(s*N), nrow = N, ncol = s)
     volterra <- t(apply(dW, 1, function(row) L%*%row))
     
     add0(volterra)
}

sim_volterra_eigen    <- function(H, n, TT = 1, N = 1) {
     s     <- ceiling(n*TT)
     Sigma <- covmat_volterra(H, n, TT)
     fbm   <- mvrnorm(n = N, mu = rep(0,s), Sigma = Sigma)
     add0(volterra)
}

sim_volterra <- function(H, n, TT = 1, N = 1) sim_volterra_cholesky(H, n, TT, N) 

sim_bm_volterra_cholesky <- function(H, rho, n, TT = 1, N = 1) {
     s     <- ceiling(n*TT)
     
     Sigma <- covmat_bm_volterra(H, rho, n, TT)
     L     <- t(chol(Sigma))
     dW    <- matrix(rnorm(2*s*N), nrow = N, ncol = 2*s)
     joint <- t(apply(dW, 1, function(row) L%*%row))
     
     if (N == 1) bm <- add0(joint[1:s])
     else        bm <- add0(joint[,1:s]) 
     
     if (N == 1) volterra <- add0(joint[(s+1):(2*s)])
     else        volterra <- add0(joint[,(s+1):(2*s)]) 
     
     return(list(bm = bm, volterra = volterra))
}

sim_bm_volterra_eigen <- function(H, rho, n, TT = 1, N = 1) {
     s     <- ceiling(n*TT)
     Sigma <- covmat_bm_volterra(H, rho, n, TT)
     joint   <- mvrnorm(n = N, mu = rep(0,2*s), Sigma = Sigma)
     
     if (N == 1) bm <- add0(joint[1:s])
     else        bm <- add0(joint[,1:s]) 
     
     if (N == 1) volterra <- add0(joint[(s+1):(2*s)])
     else        volterra <- add0(joint[,(s+1):(2*s)]) 
     
     return(list(bm = bm, volterra = volterra))
}

sim_bm_volterra <- function(H, rho, n, TT = 1, N = 1) sim_bm_volterra_cholesky(H, rho, n, TT, N) 
