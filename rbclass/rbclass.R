###########################
### IMPORTS & LIBRARIES ###
###########################

source("generalclass/simclass.R")
source("rbclass/rbsim.R")
source("rbclass/rbpricing.R")
source("rbclass/rbcalibrate.R")

#############
### CLASS ###
#############

roughbergomiclass <- function(n = 100, N = 1000, TT = NULL, a = -0.43, rho = -0.9, eta = 1.9, xi = 0.235^2, seed = 0) {
     
     rbclass        <- simulateclass(n = n, N = N, TT = TT, seed = seed)
     class(rbclass) <- c("roughbergomiclass","simulateclass")
     
     rbclass$paths     <- pathsclass()
     rbclass$vars      <- varclass(c("a","rho","eta","xi"), c(a, rho, eta, xi))
     rbclass$varbounds <- list(lb = c(-0.49,-0.99,0.001,0.001), ub = c(0.49,0.99,5,5))

     return(rbclass)
}

############################
##### FUNCTIONS TO USE #####
############################

simulate.roughbergomiclass  <- function(...) simulate_rb_antimixed(...)
price.roughbergomiclass     <- function(...) price_rb_mixed(...)
calibrate.roughbergomiclass <- function(...) calibraterb_mixed(...)

