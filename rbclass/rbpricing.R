###########################
### IMPORTS & LIBRARIES ###
###########################

########################
##### TURBOCHARGED #####
########################

price_rb_mixed <- function(sim, estimator = c("mixed", "conditional", "controlled")) {
     sim$priceinfo$starttime <- Sys.time()
     
     ##### Local variables #####
     K    <- exp(sim$simgrid$k)
     ot   <- sim$simgrid$t # option price times
     steps<- sim$simgrid$s # t[steps] = ot is option price times
     dt   <- sim$timegrid$dt
     rho  <- sim$vars$rho
     
     ##### Mixed estimator #####
     mixed_vols   <- matrix(data = NA, nrow = length(K), ncol = length(steps))
     mixed_prices <- matrix(data = NA, nrow = length(K), ncol = length(steps))

     QV <- as.matrix( t( apply(sim$paths$V, 1, cumsum) * dt )[,steps] )  # [k,t]
     Q  <- apply(QV, 2, max) + 1e-9                                      # [t]

     for (j in 1:length(steps)) {
          
          if (match.arg(estimator) == "mixed" && rho != 0) {
               X    <- as.matrix(sapply(K, function(x) vec_bs(sim$paths$S1[,steps[j]], x, (1-rho^2)*QV[,j])))      # [N,k]
               Y    <- as.matrix(sapply(K, function(x) vec_bs(sim$paths$S1[,steps[j]], x, rho^2*(Q[j] - QV[,j])))) # [N,k]
               eY   <- vec_bs(1, K, rho^2*Q[j]) # [k]
               c    <- sapply(1:length(K), function(i) -cov(X[,i],Y[,i]) / var(Y[,i])) # [k]
               
               # Payoffs, prices and volatilities
               if (length(K) == 1) {
                    mixed_payoff     <- X + c*(Y-eY) # [N,1]
                    mixed_prices[,j] <- mean(mixed_payoff) # [1]
                    mixed_vols[,j]   <- bsinv(P = mixed_prices[,j], Fwd = 1, K = K, TT = ot[j]) # [1]
               } else {
                    diff             <- t(apply(Y, 1, function(x) x-eY))  # [N,k]
                    mixed_payoffs    <- X + t(apply(diff, 1, function(x) c*x)) # [N,k]
                    mixed_prices[,j] <- apply(mixed_payoffs, 2, mean) # [k]
                    mixed_vols[,j]   <- vec_bsinv(P = mixed_prices[,j], Fwd = 1, K = K, TT = ot[j]) # [k]
               }
          } else if (match.arg(estimator) == "conditional" || rho == 0) {
               X    <- as.matrix(sapply(K, function(x) vec_bs(sim$paths$S1[,steps[j]], x, (1-rho^2)*QV[,j])))      # [N,k]
               if (length(K) == 1) {
                    mixed_prices[,j] <- mean(X) # [1]
                    mixed_vols[,j]   <- bsinv(P = mixed_prices[,j], Fwd = 1, K = K, TT = ot[j]) # [1]
               } else {
                    mixed_prices[,j] <- apply(X, 2, mean) # [k]
                    mixed_vols[,j]   <- vec_bsinv(P = mixed_prices[,j], Fwd = 1, K = K, TT = ot[j]) # [k]
               }
          } else if (match.arg(estimator) == "controlled" || rho == 1) {
               X    <- as.matrix(sapply(K, function(x) sim$paths$S[,steps[j]] - x ))      # [N,k]
               X    <- pmax(X,0)
               Y    <- as.matrix(sapply(K, function(x) vec_bs(sim$paths$S[,steps[j]], x, Q[j] - QV[,j]))) # [N,k]
               eY   <- vec_bs(1, K, Q[j]) # [k]
               c    <- sapply(1:length(K), function(i) -cov(X[,i],Y[,i]) / var(Y[,i])) # [k]
               
               # Payoffs, prices and volatilities
               if (length(K) == 1) {
                    mixed_payoff     <- X + c*(Y-eY) # [N,1]
                    mixed_prices[,j] <- mean(mixed_payoff) # [1]
                    mixed_vols[,j]   <- bsinv(P = mixed_prices[,j], Fwd = 1, K = K, TT = ot[j]) # [1]
               } else {
                    diff             <- t(apply(Y, 1, function(x) x-eY))  # [N,k]
                    mixed_payoffs    <- X + t(apply(diff, 1, function(x) c*x)) # [N,k]
                    mixed_prices[,j] <- apply(mixed_payoffs, 2, mean) # [k]
                    mixed_vols[,j]   <- vec_bsinv(P = mixed_prices[,j], Fwd = 1, K = K, TT = ot[j]) # [k]
               }
          }
     }
     
     sim$simgrid$impvol   <- mixed_vols
     sim$simgrid$prices   <- mixed_prices

     sim$priceinfo$endtime <- Sys.time()
     return(sim)
}
