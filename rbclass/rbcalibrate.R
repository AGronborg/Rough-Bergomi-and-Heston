###########################
### IMPORTS & LIBRARIES ###
###########################

library(stats) # for optim

###################################
##### CALIBRATE RBCLASS MIXED #####
###################################
changedvars <- function(vclass, lastvars) {
     if ( !any(as.numeric(vclass) != as.numeric(lastvars))) return("") 
     names(lastvars)[as.numeric(vclass) != as.numeric(lastvars)]
}

calibrate_rb_mixed <- function(simclass, plottrace = FALSE, trackpars = FALSE, minimize = c("impvol","prices"), control = list(), simfunc = simulate, pricefunc = price, ...) {
     simclass$calinfo$starttime <- Sys.time()
     
    opfunc <- function(vars) {
          changed <- changedvars(vars,lastvars)
          if (!("a" %in% changed) && !("rho" %in% changed)) {
            simclass <<- setvars(simclass, names = names(calvars), values = vars)
            simclass <<- simulate_rb_antimixed(simclass, skip = c("W1","Y","W2","B","S"))
          } else {
            simclass <<- setvars(simclass, names = names(calvars), values = vars)
            simclass <<- simulate_rb_antimixed(simclass)
          }
          lastvars <<- setvars.varclass(lastvars, names = names(lastvars), values = vars)
          simclass <<- price_rb_mixed(simclass)
          if (is.null(simclass$simgrid$prices) && minimize == "prices") simclass$simgrid$prices <- getprices(simclass$simgrid)

          sim <- getsmiles(simclass$simgrid, simclass$empgrid$t, simclass$empgrid$k, minimize)
          emp <- getsmiles(simclass$empgrid, type = minimize)

          se <- sum((sim-emp)^2*simclass$calweights)
          
          if (plottrace) plot(simclass, type = minimize)
          if (trackpars) print(c(vars,se))
          
          return(se)
    }
    
     calvars  <- identifyvars(names(simclass$vars), ...) # vars to calibrate
     lastvars <- setvars(calvars, names(calvars), rep(-1000,length(calvars)))
     bounds   <- getbounds(simclass, names(calvars))
     simclass <- simulate_rb_antimixed(simclass)
     
     minimize <- match.arg(minimize)
     if (minimize == "prices" && is.null(simclass$empgrid$prices)) simclass$empgrid$prices <- getprices(simclass$empgrid)
     
     op <- optim(par = as.numeric(calvars), fn = opfunc, lower = bounds$lb, upper = bounds$ub, method = "L-BFGS-B", control = control)
     
     simclass$calinfo$op      <- op
     simclass$calinfo$calvars <- names(calvars)
     
     simclass <- setvars(simclass, names = names(calvars), values = op$par)
     simclass <- simulate_rb_antimixed(simclass)
     simclass <- price_rb_mixed(simclass)
     
     simclass$calinfo$endtime <- Sys.time()
     return(simclass)
}

