###########################
### IMPORTS & LIBRARIES ###
###########################

library(stats)  # uniroot
library(MASS)   # mvrnorm

####################
##### SIMULATE #####
####################

simulate_rb <- function(rbclass, skip = "", antithetic = FALSE) {
     rbclass$siminfo$starttime <- Sys.time()
     rbclass <- setseed(rbclass)
     
     N   <- rbclass$N            # paths
     if (antithetic) N <- N/2    # N must be even if antithetic
     
     n   <- rbclass$timegrid$n   # time steps per year 
     s   <- rbclass$timegrid$s   # time steps
     dt  <- rbclass$timegrid$dt
     t   <- rbclass$timegrid$t   # vector of times
     
     a   <- rbclass$vars$a
     xi  <- rbclass$vars$xi
     eta <- rbclass$vars$eta
     rho <- rbclass$vars$rho
     e   <- c(0,0)
     c   <- cov_hybrid(a,n)
     
     if (!("W1" %in% skip)) {
          dW1 <- fdW1(e, c, N, s)
          rbclass$paths$dW11 <- dW1[,,1]
          rbclass$paths$dW12 <- dW1[,,2]
     }
     if (!("Y"  %in% skip)) rbclass$paths$Y   <- fY(rbclass$paths$dW11, rbclass$paths$dW12, N, s, a, n) 
     if (!("V"  %in% skip)) {
          if (antithetic) {
               rbclass$paths$V               <- matrix(NA, nrow = 2*N, ncol = s + 1)
               rbclass$paths$V[1:N,]         <- fV( rbclass$paths$Y, t, a, xi, eta)
               rbclass$paths$V[(N+1):(2*N),] <- fV(-rbclass$paths$Y, t, a, xi, eta)
          } else {
               rbclass$paths$V  <- fV(rbclass$paths$Y, t, a, xi, eta)
          }
     }
     if (!("S1" %in% skip)) {
          if (antithetic) {
               rbclass$paths$S1               <- matrix(NA, nrow = 2*N, ncol = s + 1)
               rbclass$paths$S1[1:N,]         <- fS1(rbclass$paths$V[1:N,],          rbclass$paths$dW11, rho, dt) 
               rbclass$paths$S1[(N+1):(2*N),] <- fS1(rbclass$paths$V[(N+1):(2*N),], -rbclass$paths$dW11, rho, dt) 
          } else {
               rbclass$paths$S1  <- fS1(rbclass$paths$V, rbclass$paths$dW11, rho, dt) 
          }
     }
     
     if (!("W2" %in% skip)) rbclass$paths$dW2 <- fdW2(N, s, dt)
     if (!("B"  %in% skip)) rbclass$paths$dB  <- fdB(rbclass$paths$dW11, rbclass$paths$dW2, rho)
     if (!("S"  %in% skip)) rbclass$paths$S   <- fS(rbclass$paths$V, rbclass$paths$dB, dt) 
     
     rbclass$siminfo$endtime <- Sys.time()
     return(rbclass)
}

simulate_rb_slow <- function(rbclass, skip = "", antithetic = FALSE) { 
     # slower than simulate_rb when doing antithetics, but give all paths
     rbclass$siminfo$starttime <- Sys.time()
     rbclass <- setseed(rbclass)
     
     N   <- rbclass$N            # paths
     if (antithetic) N <- N/2    # N must be even if antithetic
     if (antithetic) skip <- c("W2","B","S") # Antithetic only works with mixed
     
     n   <- rbclass$timegrid$n   # time steps per year 
     s   <- rbclass$timegrid$s   # time steps
     dt  <- rbclass$timegrid$dt
     t   <- rbclass$timegrid$t   # vector of times
     
     a   <- rbclass$vars$a
     xi  <- rbclass$vars$xi
     eta <- rbclass$vars$eta
     rho <- rbclass$vars$rho
     e   <- c(0,0)
     c   <- cov_hybrid(a,n)
     
     if (!("W1" %in% skip)) {
          dW1 <- fdW1(e, c, N, s)
          if (antithetic) {
               rbclass$paths$dW11 <- rbind(dW1[,,1],-dW1[,,1])
               rbclass$paths$dW12 <- rbind(dW1[,,2],-dW1[,,2])
          } else {
               rbclass$paths$dW11 <- dW1[,,1]
               rbclass$paths$dW12 <- dW1[,,2]
          }
     }
     if (!("Y"  %in% skip))  {
          if (antithetic) {
               Y <- fY(rbclass$paths$dW11[1:N,], rbclass$paths$dW12[1:N,], N, s, a, n)
               rbclass$paths$Y <- rbind(Y,-Y)
          } else {
               rbclass$paths$Y <- fY(rbclass$paths$dW11, rbclass$paths$dW12, N, s, a, n)
          }
     }
     if (!("V"  %in% skip)) rbclass$paths$V  <- fV(rbclass$paths$Y, t, a, xi, eta)
     if (!("S1" %in% skip)) rbclass$paths$S1 <- fS1(rbclass$paths$V, rbclass$paths$dW11, rho, dt) 
     
     if (!("W2" %in% skip)) rbclass$paths$dW2 <- fdW2(N, s, dt)
     if (!("B"  %in% skip)) rbclass$paths$dB  <- fdB(rbclass$paths$dW11, rbclass$paths$dW2, rho)
     if (!("S"  %in% skip)) rbclass$paths$S   <- fS(rbclass$paths$V, rbclass$paths$dB, dt) 

     rbclass$siminfo$endtime <- Sys.time()
     return(rbclass)
}
          
simulate_rb_standard            <- function(rbclass, skip = "") simulate_rb(rbclass, skip = skip, antithetic = FALSE)
simulate_rb_mixed               <- function(rbclass, skip = c("W2","B","S")) return(simulate_rb(rbclass, skip = skip, antithetic = FALSE))
simulate_rb_antimixed           <- function(rbclass, skip = c("W2","B","S")) return(simulate_rb(rbclass, skip = skip, antithetic = TRUE))
simulate_rb_antimixed_withpaths <- function(rbclass, skip = c("W2","B","S")) return(simulate_rb_slow(rbclass, skip = skip, antithetic = TRUE))

#################
##### UTILS #####
#################

##### Kernel #####
g <- function(x,a) (x^a)

##### Discretization (minimising hybrid scheme error) #####
b <- function(k,a) ((k^(a+1)-(k-1)^(a+1))/(a+1))^(1/a)

##### Covariance matrix for hybrid scheme (assuming kappa = 1) #####
cov_hybrid <- function(a, n) {
     cov <- matrix(rep(0,4), 2, 2)
     cov[1,1] <- 1/n
     cov[2,1] <- 1/((a+1)*n^(a+1))
     cov[1,2] <- cov[2,1]
     cov[2,2] <- 1/((2*a+1)*n^(2*a+1))
     
     return (cov)
}

#################
##### PATHS #####
#################

##### dW1 #####
fdW1 <- function(e, c, N, s) {
     paths <- array(data = NA, dim = c(N,s,2)) # (N,s,W)
     
     for (i in 1:N) {
          paths[i,,]  <- mvrnorm(n = s, mu = e, Sigma = c)
     }
     
     return(paths)
}

##### Y #####
fY <- function(dW11, dW12, N, s, a, n) {
     
     # Volterra process
     Y1 <- matrix(data = 0, nrow = N, ncol = s + 1) # Exact integrals
     Y2 <- matrix(data = 0, nrow = N, ncol = s + 1) # Riemann sums
     
     # Exact integral
     for (i in 2:(s+1)) Y1[,i] <- dW12[,i-1] # Assumes kappa = 1
     
     # Arrays for convolution
     G <- rep(0,s + 1) # Gamma
     
     for (k in 2:s) G[k+1] <- g( b(k,a)/n, a )
     
     X <- dW11[,] # Xi
     if (!is.matrix(X)) X <- as.matrix(X)
     
     GX <- matrix(data = 0, nrow = N, ncol = length(X[1,])+length(G)-1) # row = paths, col = 2 * timesteps
     
     for (i in 1:N) {
          GX[i,] <- convolve(G, rev(X[i,]), type="o") # convolve(x,rev(y),type="o") is normal convolution of x and y
     }
     
     Y2 <- GX[,1:(1+s)]
     
     Y <- sqrt(2*a + 1) * (Y1+Y2)
     
     return(Y)
}

##### dW2 #####
fdW2 <- function(N,s,dt) {
     dW <- matrix(data = rnorm(n = N*s, sd = sqrt(dt)), nrow = N, ncol = s)
     return(dW)
}

##### dB #####
fdB <- function(dW11, dW2, rho) {
     dB <- rho * dW11[,] + sqrt(1 - rho^2) * dW2
     return (dB)
}

##### V #####
fV <- function(Y, t, a, xi, eta) {
     m    <- eta*Y
     v    <- 0.5*eta^2*t^(2*a + 1)
     diff <- apply(m, 1, function(x) x - v)
     V    <- t( xi*exp(diff) )
     return(V)
}

##### S #####
fS <- function(V, dB, dt, S0 = 1) {
     increments <- sqrt(V[,-ncol(V)])*dB - 0.5*V[,-ncol(V)]*dt
     if (!is.matrix(increments)) increments <- as.matrix(increments)
     
     integral <- t( apply(increments, 1, cumsum) )
     
     S <- matrix(data = 0, nrow = nrow(V), ncol = ncol(V))
     
     S[,1] <- S0
     S[,2:ncol(V)] <- S0 * exp(integral)
     
     return(S)
}

##### S1 #####
fS1 <- function(V, dW11, rho, dt, S0 = 1) {
     increments <- rho * sqrt(V[,-ncol(V)]) * dW11[,] - 0.5 * rho^2 * V[,-ncol(V)] * dt
     if (!is.matrix(increments)) increments <- as.matrix(increments)
     
     integral <- t( apply(increments, 1, cumsum) )
     
     S <- matrix(data = 0, nrow = nrow(V), ncol = ncol(V))
     
     S[,1] <- S0
     S[,2:ncol(V)] <- S0 * exp(integral)
     
     return(S)
}

##########################################
##### ALTERNATIVE ANTITHETIC VERSION #####
##########################################
# 
# simulate_rb <- function(rbclass, skip = "", antithetic = FALSE) {
#      
#      rbclass <- setseed(rbclass)
#      
#      N   <- rbclass$N            # paths
#      if (antithetic) N <- N/2    # N must be even if antithetic
#      if (antithetic) skip <- c("W2","B","S") # Antithetic only works with mixed
#      
#      n   <- rbclass$timegrid$n   # time steps per year 
#      s   <- rbclass$timegrid$s   # time steps
#      dt  <- rbclass$timegrid$dt
#      t   <- rbclass$timegrid$t   # vector of times
#      
#      a   <- rbclass$vars$a
#      xi  <- rbclass$vars$xi
#      eta <- rbclass$vars$eta
#      rho <- rbclass$vars$rho
#      e   <- c(0,0)
#      c   <- cov_hybrid(a,n)
#      
#      rbclass$siminfo$starttime <- Sys.time()
#      if (!("W1" %in% skip)) {
#           dW1 <- fdW1(e, c, N, s)
#           if (antithetic) {
#                rbclass$paths$dW11 <- rbind(dW1[,,1],-dW1[,,1])
#                rbclass$paths$dW12 <- rbind(dW1[,,2],-dW1[,,2])
#           } else {
#                rbclass$paths$dW11 <- dW1[,,1]
#                rbclass$paths$dW12 <- dW1[,,2]
#           }
#      }
#      if (!("Y"  %in% skip))  {
#           if (antithetic) {
#                Y <- fY(rbclass$paths$dW11[1:N,], rbclass$paths$dW12[1:N,], N, s, a, n)
#                rbclass$paths$Y <- rbind(Y,-Y)
#           } else {
#                rbclass$paths$Y <- fY(rbclass$paths$dW11, rbclass$paths$dW12, N, s, a, n)
#           }
#      }
#      if (!("V"  %in% skip)) rbclass$paths$V  <- fV(rbclass$paths$Y, t, a, xi, eta)
#      if (!("S1" %in% skip)) rbclass$paths$S1 <- fS1(rbclass$paths$V, rbclass$paths$dW11, rho, dt) 
#      
#      if (!("W2" %in% skip)) rbclass$paths$dW2 <- fdW2(N, s, dt)
#      if (!("B"  %in% skip)) rbclass$paths$dB  <- fdB(rbclass$paths$dW11, rbclass$paths$dW2, rho)
#      if (!("S"  %in% skip)) rbclass$paths$S   <- fS(rbclass$paths$V, rbclass$paths$dB, dt) 
#      rbclass$siminfo$endtime <- Sys.time()
#      
#      return(rbclass)
# }