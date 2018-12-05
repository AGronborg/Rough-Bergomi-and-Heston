###########################
### IMPORTS & LIBRARIES ###
###########################

library(stats)  # uniroot
library(MASS)   # mvrnorm

####################
##### SIMULATE #####
####################

simulate_rb <- function(rbclass, skip = "", antithetic = FALSE, exact = FALSE) {
     rbclass$siminfo$starttime <- Sys.time()
     rbclass <- setseed(rbclass)
     
     N   <- rbclass$N            # paths
     if (antithetic) N <- N/2    # N must be even if antithetic
     if (!(N %% 1 == 0)) stop("N has to be even when using antithetics")

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
     if (!("Y"  %in% skip)) rbclass$paths$Y <- fY(rbclass$paths$dW11, rbclass$paths$dW12, N, s, a, n)
     if (!("V"  %in% skip)) {
          rbclass$paths$V  <- fV(rbclass$paths$Y, t, a, xi, eta)
          if (antithetic) rbclass$paths$V <- rbind(rbclass$paths$V, fV(-rbclass$paths$Y, t, a, xi, eta))
     }
     if (!("S1" %in% skip)) {
          rbclass$paths$S1  <- fS1(rbclass$paths$V[1:N,], rbclass$paths$dW11, rho, dt)
          if (antithetic) rbclass$paths$S1 <- rbind(rbclass$paths$S1, fS1(rbclass$paths$V[(N+1):(2*N),], -rbclass$paths$dW11, rho, dt) )
     }
     
     if (!("W2" %in% skip)) rbclass$paths$dW2 <- fdW2(N, s, dt)
     if (!("B"  %in% skip)) rbclass$paths$dB  <- fdB(rbclass$paths$dW11, rbclass$paths$dW2, rho)
     if (exact) {
          joint            <- sim_bm_volterra(a + 0.5, rho, n, s/n, N)
          rbclass$paths$Y  <- joint$volterra
          rbclass$paths$dB <- t(diff(t(joint$bm)))
          rbclass$paths$V  <- fV(rbclass$paths$Y, t, a, xi, eta)
     }
     if (!("S"  %in% skip)) rbclass$paths$S   <- fS(rbclass$paths$V, rbclass$paths$dB, dt) 
     
     rbclass$siminfo$endtime <- Sys.time()
     return(rbclass)
}

simulate_rb_exact     <- function(rbclass, skip = c("W1", "W2", "B", "Y", "V", "S1")) simulate_rb(rbclass, skip = skip, antithetic = FALSE, exact = TRUE)
simulate_rb_standard  <- function(rbclass, skip = "S1") simulate_rb(rbclass, skip = skip, antithetic = FALSE, exact = FALSE)
simulate_rb_mixed     <- function(rbclass, skip = c("W2","B","S")) simulate_rb(rbclass, skip = skip, antithetic = FALSE, exact = FALSE)
simulate_rb_antimixed <- function(rbclass, skip = c("W2","B","S")) simulate_rb(rbclass, skip = skip, antithetic = TRUE, exact = FALSE)

simulate_rb_antimixed_withpaths <- function(rbclass, skip = c("W2","B","S")) {
     rbclass <- simulate_rb_antimixed(rbclass, skip = skip)
     rbclass <- rb_add_antithetic_paths(rbclass)
     rbclass$siminfo$endtime <- Sys.time()
     return(rbclass)
}

rb_add_antithetic_paths <- function(rbclass) {
     rbclass$paths$dW11 <- rbind(rbclass$paths$dW11, -rbclass$paths$dW11)
     rbclass$paths$dW12 <- rbind(rbclass$paths$dW12, -rbclass$paths$dW12)
     rbclass$paths$Y    <- rbind(rbclass$paths$Y,    -rbclass$paths$Y)
     return(rbclass)
}

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
fdW1 <- function(mu, Sigma, N, s) {
     dW <- array(data = NA, dim = c(N,s,length(mu))) # (N,s,W)
     
     for (i in 1:N) dW[i,,]  <- mvrnorm(n = s, mu = mu, Sigma = Sigma)
     
     return(dW)
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