###########################
### IMPORTS & LIBRARIES ###
###########################

source("generalclass/utils.R")

#################
##### PRICE #####
#################

price               <- function(simclass, ...) UseMethod("price", simclass)
price.simulateclass <- function(...) price_standard(...)

##########################
##### PRICE SIMCLASS #####
##########################

price_standard <- function(simclass) {
     simclass$priceinfo$starttime <- Sys.time()
     ##### Local variables #####
     K    <- exp(simclass$simgrid$k) # strikes
     t    <- simclass$simgrid$t      # option price times
     s    <- simclass$simgrid$s      # t[s] is option price times
     
     ##### Base estimator #####
     prices <- matrix(data = NA, nrow = length(K), ncol = length(s))
     impvol <- matrix(data = NA, nrow = length(K), ncol = length(s))
     
     for (j in 1:length(s)) {
          Sj             <- simclass$paths$S[,s[j]]
          diff           <- t( sapply(Sj, function(x) x-K) )
          payoffs        <- pmax(diff, 0)
          
          if (length(t) == 1) prices[,j] <- mean(payoffs)
          else if (length(t) > 1) prices[,j] <- apply(payoffs, 2, mean)
          impvol[,j]     <- vec_bsinv(P = prices[,j], Fwd = 1, K = K, TT = t[j])
     }
     
     simclass$simgrid$prices   <- prices
     simclass$simgrid$impvol   <- impvol
     
     simclass$priceinfo$endtime <- Sys.time()
     return(simclass)
}

single_price <- function(simclass, t, k, n = NULL, N = NULL, varnames = NULL, values = NULL, simfunc = simulate, pricefunc = price, type = c("impvol","price"), seed = -1, ...) {
     if (is.null(n) && !is.null(timegrid)) n <- simclass$n
     else if (is.null(n)) n <- simclass$timegrid$n
     
     simclass <- setvars(simclass, names = varnames, values = values, N = N, seed = seed, ...)
     simclass <- changetimegrid(simclass, TT = t, n = n, reset = TRUE)
     simclass <- setsimgrid(simclass, t = t, k = k)
     simclass <- simfunc(simclass)
     simclass <- pricefunc(simclass)
     if (match.arg(type) == "impvol") return(simclass$simgrid$impvol[1,1])
     else if (match.arg(type) == "price") {
          if (!is.null(simclass$simgrid$prices)) return(simclass$simgrid$prices[1,1])
          else return(bs(Fwd = 1, K = exp(k), V = ceiling(n*t)/n*simclass$simgrid$impvol[1,1]^2))
     }
}

single_price_vec <- function(simclass, t, k, n = NULL, N = NULL, ...) {
     
     if (is.null(n)) n <- simclass$timegrid$n
     if (is.null(N)) N <- simclass$N
     
     prevecfunc <- function(t, k, n, N) single_price(simclass = simclass, t = t, k = k, n = n, N = N, ...)
     vecfunc <- Vectorize(prevecfunc)
     vecfunc(t, k, n, N)
}

pricetime <- function(object, ...) UseMethod("pricetime", object)
pricetime.simulateclass <- function(simclass, units = "auto", digits = 0) gettime(simclass$priceinfo, units, digits)

simpricetime <- function(object, ...) UseMethod("simpricetime")
simpricetime.simulateclass <- function(simclass, units = "auto", digits = 0) {
     gettime(simclass$siminfo, units, digits) + gettime(simclass$priceinfo, units, digits)
}