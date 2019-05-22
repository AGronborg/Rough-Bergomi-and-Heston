###########################
### IMPORTS & LIBRARIES ###
###########################

library(MASS)   # mvrnorm

################################
##### SIMULATE HESTONCLASS #####
################################

simulate_heston <- function(hclass, scheme = 2) {
     
     hclass$siminfo$starttime <- Sys.time()
     hclass <- setseed(hclass)
     
     N  <- hclass$N
     s  <- hclass$timegrid$s
     dt <- hclass$timegrid$dt
     
     lambda <- hclass$vars$lambda
     vbar   <- hclass$vars$vbar
     v0     <- hclass$vars$v0
     eta    <- hclass$vars$eta
     rho    <- hclass$vars$rho
     
     dW <- hfdW_rhocor(dt,rho,N,s)
     V  <- hfV(dW[,,2], dt, lambda, vbar, v0, eta, scheme)
     S  <- hfS(dW[,,1], V, dt)
     

     hclass$paths$dW1 <- dW[,,1]
     hclass$paths$dW2 <- dW[,,2]
     hclass$paths$V   <- V
     hclass$paths$S   <- S
     
     hclass$siminfo$endtime <- Sys.time()
     return(hclass)
}

#################
##### PATHS #####
#################

##### dW1 and dW2 #####
hfdW <- function(mu, Sigma, N, s) {
     dW <- array(data = NA, dim = c(N,s,length(mu))) # (N,s,W)
     
     for (i in 1:N) dW[i,,]  <- mvrnorm(n = s, mu = mu, Sigma = Sigma)

     return(dW)
}
hfdW_rhocor <- function(var, rho, N, s) hfdW(c(0,0), matrix(c(var,rho*var,rho*var,var),2,2),N,s)

##### V #####
hfV <- function(dW, dt, lambda, vbar, v0, eta, scheme = 2) {
     
     V <- matrix(data = v0, nrow = nrow(dW), ncol = ncol(dW) + 1)# (N,s)

     for (j in 2:(ncol(dW)+1))
     {
          if ((scheme == 3) & (2*lambda*vbar/(eta^2) <= 1)) {
               warning("Because of the choice of lambda, vbar and eta, the variance can be negative. Thus Reflection + Milstein method is used instead.")
               scheme = 2
          }
          
          if (scheme == 0){
               ## Euler Gatheral (2.17) + Absorbtion
               V[,j] <- V[,j-1] + lambda*(vbar - V[,j-1])* dt + eta * sqrt(V[,j-1]) * dW[,j-1]
               V[V[,j] < 0,j]<- 0
          }
          else if (scheme == 1){
               ## Euler Gatheral (2.17) + Reflection
               V[,j] <- V[,j-1] + lambda*(vbar - V[,j-1])* dt + eta * sqrt(V[,j-1]) * dW
               V[,j] <- ifelse(V[,j]<0, -V[,j], V[,j])
          }
          else if (scheme == 2) {
               ## Milstein Gatheral (2.18) + Reflection
               V[,j] <- (sqrt(V[,j-1]) + eta/2*dW[,j-1])^2 - lambda*(V[,j-1]-vbar)*dt - eta^2/4*dt
               V[,j] <- ifelse(V[,j]<0, -V[,j], V[,j])     
          }
          else if (scheme == 3) {
               ## Alfonsi (Gatheral p.23)
               V[,j] <- V[,j-1] -lambda*(v-vbar)*dt +eta*sqrt(V[,j-1])*dW - eta^2/2*dt      
          }
     }
     
     return(V)
}

##### S #####
hfS <- function(dW,V,dt) {
     
     S <- matrix(data = 1, nrow = nrow(V), ncol = ncol(V))# (N,s)
     
     for (j in 2:(ncol(dW)+1)) {
          S[,j] <- S[,j-1]*exp(-V[,j-1]/2*dt + sqrt(V[,j-1]) * dW[,j-1])
     }
     return(S)
}