###########################
### IMPORTS & LIBRARIES ###
###########################

library(stats)  # uniroot
library(MASS)   # mvrnorm
library(hypergeo) # hypergeo

####################
##### SIMULATE #####
####################

simulate_rb <- function(rbclass, skip = "", antithetic = FALSE, exact = FALSE, kappa = 1) {
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
     
     if (!("W1" %in% skip)) {
          e   <- rep(0,kappa+1)
          c   <- cov_hybrid(a,n,kappa)
          dW1 <- fdW1(e, c, N, s)
          
          rbclass$paths$dW11 <- dW1[,,1]
          rbclass$paths$dW12 <- dW1[,,-1,drop=FALSE]
     }
     if (!("Y"  %in% skip)) rbclass$paths$Y <- fY(rbclass$paths$dW11, rbclass$paths$dW12, N, s, a, n, kappa)
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

simulate_rb_exact     <- function(rbclass, skip = c("W1", "W2", "B", "Y", "V", "S1"), ...) simulate_rb(rbclass, skip = skip, antithetic = FALSE, exact = TRUE, ...)
simulate_rb_standard  <- function(rbclass, skip = "S1", ...) simulate_rb(rbclass, skip = skip, antithetic = FALSE, exact = FALSE, ...)
simulate_rb_mixed     <- function(rbclass, skip = c("W2","B","S"), ...) simulate_rb(rbclass, skip = skip, antithetic = FALSE, exact = FALSE, ...)
simulate_rb_antimixed <- function(rbclass, skip = c("W2","B","S"), ...) simulate_rb(rbclass, skip = skip, antithetic = TRUE, exact = FALSE, ...)

simulate_rb_antimixed_withpaths <- function(rbclass, skip = c("W2","B","S"), ...) {
     rbclass <- simulate_rb_antimixed(rbclass, skip = skip, ...)
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
g_kernel <- function(x,a) (x^a)
g_kernel <- Vectorize(g_kernel)

##### Discretization (minimising hybrid scheme error) #####
b <- function(k,a) ((k^(a+1)-(k-1)^(a+1))/(a+1))^(1/a)
b_weights <- function(k,a) ((k^(a+1)-(k-1)^(a+1))/(a+1))^(1/a)
b_weights <- Vectorize(b_weights)

##### Covariance matrix for hybrid scheme #####
cov_hybrid <- function(a, n, kappa = 1) {
     cov <- matrix(NA, kappa + 1, kappa + 1)
     cov[1,1] <- 1/n
     
     for (j in 2:(kappa + 1))
          cov[1,j] <- cov[j,1] <- ((j-1)^(a+1)-(j-2)^(a+1))/((a+1)*n^(a+1))
     for (j in 2:(kappa + 1))
          cov[j,j] <- ((j-1)^(2*a+1)-(j-2)^(2*a+1))/((2*a+1)*n^(2*a+1))
     
     if (kappa > 1) {
          for (j in 2:kappa) {
               for (k in (j+1):(kappa + 1)) {
                    term1 <- (j-1)^(a+1)*(k-1)^a*hypergeo(-a,1,a+2,((j-1)/(k-1)))
                    term2 <- (j-2)^(a+1)*(k-2)^a*hypergeo(-a,1,a+2,((j-2)/(k-2)) )
                    cov[j,k] <- cov[k,j] <- 1/((a+1)*n^(2*a+1))*(term1-term2)
               }
          }
     }
     
     return(Re(cov))
}

cov_hybrid(-0.43, 10, 2)

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
fY <- function(dW11, dW12, N, s, a, n, kappa = 1) {

     Y <- NULL
     if (N == 1) {
          Y1 <- c(0, dW12[,,1])
          if (kappa > 1) for (k in 2:kappa) Y1 <- Y1 + c(rep(0,k),dW12[,-((s+2-k):s),k])
     } else {
          Y1 <- cbind(0,dW12[,,1]) # Exact integrals
          if (kappa > 1) for (k in 2:kappa) Y1 <- Y1 + cbind(matrix(0, nrow = N, ncol = k),dW12[,-((s+2-k):s),k])
     }
     
     # Arrays for convolution
     Gamma <- rep(0,s + 1) # Gamma
     for (k in (kappa+1):s) Gamma[k+1] <- g_kernel( b_weights(k,a)/n, a )
     
     Xi <- dW11 # Xi
     if (N == 1) Xi <- t(as.matrix(Xi))
     
     GX <- matrix(data = 0, nrow = N, ncol = length(Xi[1,])+length(Gamma)-1) # row = paths, col = 2 * timesteps
     
     for (i in 1:N) GX[i,] <- convolve(Gamma, rev(Xi[i,]), type="o") # convolve(x,rev(y),type="o") is normal convolution of x and y
     
     Y2 <- GX[,1:(1+s)] # Riemann sums
     
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

#################
##### EXACT #####
#################

G <- function(x, H) {
     gamma <- 0.5 - H
     Re( (1-2*gamma)/(1-gamma) * x^(-gamma) * hypergeo(1, gamma, 2 - gamma, 1/x) )
}

cov_bm_volterra <- function(bm_t, volterra_t, H, rho) {
     DH <- sqrt(2*H)/(H+0.5)
     rho*DH*( volterra_t^(H+0.5) - (volterra_t-min(volterra_t,bm_t))^(H+0.5) )
}

cov_volterra <- function(s, t, H) min(s,t)^(2*H)*G(max(t/s,s/t), H)

cov_bm <- function(s, t) min(s,t)

covmat_bm_volterra <- function(H, rho, n, TT = 1) {
     s  <- ceiling(n*TT)      # time steps total
     TT <- s/n                # maturity
     dt <- 1/n                # time step
     t  <- cumsum(rep(dt, s)) # time step sequence
     
     Sigma <- matrix(data = NA, nrow = 2*s, ncol = 2*s)
     
     for (i in 1:s) for (j in i:s) 
          Sigma[i,j] <- Sigma[j,i] <- cov_bm(t[i], t[j])
     for (i in 1:s) for (j in 1:s)
          Sigma[i,s+j] <- Sigma[s+j,i] <- cov_bm_volterra(t[i], t[j], H, rho)
     for (i in 1:s) for (j in i:s) 
          Sigma[s+i,s+j] <- Sigma[s+j,s+i] <- cov_volterra(t[i], t[j], H)

     return(Sigma)
}

sim_bm_volterra <- function(H, rho, n, TT = 1, N = 1) {
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

##########################$
##### OTHER PROCESSES #####
##########################$

##### fBm #####
cov_fbm <- function(s, t, H) 0.5*(t^(2*H)+s^(2*H)-abs(t-s)^(2*H))

covmat_fbm <- function(H, n, TT = 1) {
     s  <- ceiling(n*TT)      # time steps total
     TT <- s/n                # maturity
     dt <- 1/n                # time step
     t  <- cumsum(rep(dt, s)) # time step sequence
     
     Sigma <- matrix(data = NA, nrow = s, ncol = s)
     
     for (i in 1:s) for (j in i:s)
          Sigma[i,j] <- Sigma[j,i] <- cov_fbm(t[i], t[j], H)
     
     return(Sigma)
}

sim_fbm <- function(H, n, TT = 1, N = 1) {
     s     <- ceiling(n*TT)
     
     Sigma <- covmat_fbm(H,n,TT)
     L     <- t(chol(Sigma))
     dW   <- matrix(rnorm(s*N), nrow = N, ncol = s)
     fbm  <- t(apply(dW, 1, function(row) L%*%row))
     
     add0(fbm)
}

##### VOLTERRA #####
covmat_volterra <- function(H, n, TT = 1) {
     s  <- ceiling(n*TT)      # time steps total
     TT <- s/n                # maturity
     dt <- 1/n                # time step
     t  <- cumsum(rep(dt, s)) # time step sequence
     
     Sigma <- matrix(data = NA, nrow = s, ncol = s)
     
     for (i in 1:s) for (j in i:s)
          Sigma[i,j] <- Sigma[j,i] <- cov_volterra(t[i], t[j], H)
     
     return(Sigma)
}

sim_volterra <- function(H, n, TT = 1, N = 1) {
     s     <- ceiling(n*TT)
     
     Sigma <- covmat_volterra(H, n, TT)
     L     <- t(chol(Sigma))
     dW    <- matrix(rnorm(s*N), nrow = N, ncol = s)
     volterra <- t(apply(dW, 1, function(row) L%*%row))
     
     add0(volterra)
}

sim_volterra_hybrid  <- function(H, n, TT = 1, N = 1, kappa = 1) {
     s   <- ceiling(n*TT)
     a   <- H - 0.5
     e   <- rep(0,kappa+1)
     c   <- cov_hybrid(a,n,kappa)
     
     dW1 <- fdW1(e, c, N, s)
     fY(dW1[,,1], dW1[,,-1,drop=FALSE], N, s, a, n, kappa)
}

