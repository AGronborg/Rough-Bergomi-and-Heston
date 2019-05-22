###########################
### IMPORTS & LIBRARIES ###
###########################

source("generalclass/simclass.R")
source("hestonclass/hestonsim.R")
source("hestonclass/hestonpricing.R")

#############
### CLASS ###
#############

hestonclass <- function(n = 100, N = 1000, TT = NULL, lambda = 0.2, vbar = 0.3^2, v0 = 0.2^2, eta = 0.02, rho = 0, seed = -1) {
     
     hclass <- simulateclass(n = n, N = N, TT = TT, seed = seed)
     class(hclass) <- c("hestonclass","simulateclass")

     hclass$paths      <- pathsclass()
     hclass$vars       <- varclass(c("lambda","vbar","v0","eta","rho"), c(lambda, vbar, v0, eta, rho))
     hclass$varbounds  <- list(lb = c(0.001,0.001,0.001,0.001,-0.99), ub = c(10,10,10,10,0.99))
     
     return(hclass)
}

############################
##### FUNCTIONS TO USE #####
############################

simulate.hestonclass    <- function(...) identity(...)
price.hestonclass       <- function(...) price_heston_closedform_lipton(...)
# calibrate.hestonclass <- function(...) calibrate.simulateclass(...)

