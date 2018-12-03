###########################
### IMPORTS & LIBRARIES ###
###########################

source("generalclass/utils.R")
library(fourierin)

#############################
##### PRICE HESTONCLASS #####
#############################

##### GATHERAL #####
heston_closedform_gatheral <- function(lambda, vbar, eta, rho, v0, tau, K) {
     
     integrand <- function(u, lambda, vbar, eta, rho, v0, tau, K, j) {
          x <- -log(K)
          a <- lambda * vbar
          
          b     <- lambda                - rho*eta*(j == 1)
          alpha <- -u^2/2 -1i*u/2        + 1i*u   *(j == 1)
          beta  <- lambda - 1i*rho*eta*u - rho*eta*(j == 1)
          
          gamma  <- eta^2/2
          d      <- sqrt(beta^2 - 4*alpha*gamma)
          rminus <- (beta-d)/(2*gamma)
          rplus  <- (beta+d)/(2*gamma)
          g      <- rminus/ rplus
          
          D      <- rminus * (1 - exp(-d*tau))/(1-g*exp(-d*tau))
          C      <- lambda * (rminus * tau - 2/(eta^2) * log( (1-g*exp(-d*tau))/(1-g) ) )
          
          Re(exp(C*vbar + D*v0 + 1i*u*x) / (1i * u))
     }
     
     P <- function(lambda, vbar, eta, rho, v0, tau, K, j) {
          0.5 + 1/pi * integrate(integrand, lower = 0, upper = Inf, lambda, vbar, eta, rho, v0, tau, K, j, subdivisions=1000)$value
     }
     
     P(lambda, vbar, eta, rho, v0, tau, K, 1) - K*P(lambda, vbar, eta, rho, v0, tau, K, 0)
}
heston_closedform_gatheral_vec <- Vectorize(heston_closedform_gatheral)

price_heston_closedform_gatheral <- function(hclass) {
     hclass$priceinfo$starttime <- Sys.time()
     hclass$simgrid$prices <- matrix(NA, length(hclass$simgrid$k), length(hclass$simgrid$t))
     for (j in 1:length(hclass$simgrid$t)) {
               prices <- heston_closedform_gatheral_vec(lambda   = hclass$vars$lambda,
                                             vbar     = hclass$vars$vbar,
                                             eta      = hclass$vars$eta,
                                             rho      = hclass$vars$rho,
                                             v0       = hclass$vars$v0,
                                             tau      = hclass$simgrid$t[j],
                                             K        = exp(hclass$simgrid$k))
               hclass$simgrid$prices[,j] <- prices
               hclass$simgrid$impvol[,j] <- vec_bsinv(P = prices, Fwd = 1, K = exp(hclass$simgrid$k), TT = hclass$simgrid$t[j], o = "call")
     }
     hclass$priceinfo$endtime <- Sys.time()
     return(hclass)
}

##### LIPTON #####
heston_closedform_lipton <- function(kap, theta, eps, rho, v, tau, S = 1, le = -0.5, ue = 0.5, res = 2^10, eval_grid = NULL) {
     
     integrand <- function(k) {
          kaphat    <- kap-rho*eps/2
          zeta      <- sqrt(k^2*eps^2*(1-rho^2) + 1i*2*k*eps*rho*kaphat + kaphat^2 + eps^2/4) # This is the right sqrt
          psiminus  <-  (1i*k*rho*eps+kaphat)+zeta # Really weird notation in article
          psiplus   <- -(1i*k*rho*eps+kaphat)+zeta
          beta      <- (1-exp(-zeta*tau))/(psiminus+psiplus*exp(-zeta*tau)) # Article missing -
          alpha     <- -kap*theta/eps^2 * ( psiplus*tau + 2*log((psiminus+psiplus*exp(-zeta*tau))/(2*zeta)) )
          
          1/(k^2+0.25) * exp(alpha-(k^2+0.25)*beta*v) # exp(0.5*X) to be added later
     }
     
     out <- fourierin(f = integrand, lower_int = -50, upper_int = 50, # freq-adj = -1 to have exp(-iXk)
                      lower_eval = le, upper_eval = ue,
                      const_adj = -1, freq_adj = -1, resolution = res, eval_grid = eval_grid)
     
     if (!is.null(eval_grid)) out <- list(w = eval_grid, values = out)
     as.numeric(S - S * exp(-0.5*out$w) * Re(out$values)) # Prices
}

price_heston_closedform_lipton <- function(hclass) {
     hclass$priceinfo$starttime <- Sys.time()
     hclass$simgrid$prices <- matrix(NA, length(hclass$simgrid$k), length(hclass$simgrid$t))
     for (j in 1:length(hclass$simgrid$t)) {
          prices <- heston_closedform_lipton(kap       = hclass$vars$lambda,
                                           theta     = hclass$vars$vbar,
                                           eps       = hclass$vars$eta,
                                           rho       = hclass$vars$rho,
                                           v         = hclass$vars$v0,
                                           tau       = hclass$simgrid$t[j],
                                           eval_grid = -hclass$simgrid$k)
          hclass$simgrid$prices[,j] <- prices
          hclass$simgrid$impvol[,j] <- vec_bsinv(P = prices, Fwd = 1, K = exp(hclass$simgrid$k), TT = hclass$simgrid$t[j], o = "call")
     }
     hclass$priceinfo$endtime <- Sys.time()
     return(hclass)
}
