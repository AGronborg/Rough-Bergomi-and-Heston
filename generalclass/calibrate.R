###########################
### IMPORTS & LIBRARIES ###
###########################

library(stats) # for optim

#####################
##### CALIBRATE #####
#####################

calibrate               <- function(simclass, ...) UseMethod("calibrate", simclass)
calibrate.simulateclass <- function(...) calibrate_standard(...)

##############################
##### CALIBRATE SIMCLASS #####
##############################

calibrate_standard <- function(simclass, plottrace = FALSE, trackpars = FALSE, minimize = c("impvol","prices"), control = list(), simfunc = simulate, pricefunc = price, ...) {

     simclass$calinfo$starttime <- Sys.time()
     opfunc <- function(vars) {
          simclass <<- setvars(simclass, names = names(calvars), values = vars)
          simclass <<- simfunc(simclass)
          simclass <<- pricefunc(simclass)
          if (is.null(simclass$simgrid$prices) && minimize == "prices") simclass$simgrid$prices <- getprices(simclass$simgrid)

          sim <- getsmiles(simclass$simgrid, simclass$empgrid$t, simclass$empgrid$k, minimize)
          emp <- getsmiles(simclass$empgrid, pricetype = minimize)

          se <- sum((sim-emp)^2*simclass$calweights)
          
          if (plottrace) plot(simclass, pricetype = minimize)
          if (trackpars) print(c(vars,se))
          
          return(se)
     }
     
     calvars <- identifyvars(names(simclass$vars), ...) # vars to calibrate
     bounds  <- getbounds(simclass, names(calvars))

     minimize <- match.arg(minimize)
     if (minimize == "prices" && is.null(simclass$empgrid$prices)) simclass$empgrid$prices <- getprices(simclass$empgrid)
     
     op <- optim(par = as.numeric(calvars), fn = opfunc, lower = bounds$lb, upper = bounds$ub, method = "L-BFGS-B", control = control)
     
     simclass$calinfo$op      <- op
     simclass$calinfo$calvars <- names(calvars)
     
     #simclass <- setvars(simclass, names = names(calvars), values = op$par)
     #simclass <- simfunc(simclass)
     #simclass <- pricefunc(simclass)
     
     simclass$calinfo$endtime   <- Sys.time()
     return(simclass)
}

caltime <- function(object, ...) UseMethod("caltime", object)
caltime.simulateclass <- function(simclass, units = "secs", digits = 0) gettime(simclass$calinfo, units, digits)
