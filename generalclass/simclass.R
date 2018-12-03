###########################
### IMPORTS & LIBRARIES ###
###########################

library(microbenchmark)

source("generalclass/utils.R")
source("generalclass/volgrid.R")
source("generalclass/paths.R")
source("generalclass/vars.R")
source("generalclass/simulate.R")
source("generalclass/pricing.R")
source("generalclass/calibrate.R")

###############
### CLASSES ###
###############
simulateclass <- function(n = 100, N = 1000, TT = NULL, seed = 0) {
     
     simclass <- list(timegrid     = NULL,
                      N            = N,
                      seed         = seed,
                      lastsimseed  = NULL,
                      simgrid      = NULL,
                      empgrid      = NULL,
                      calweights   = NULL,
                      calinfo      = NULL,
                      
                      vars         = NULL, # Sub classes
                      varbounds    = NULL,
                      paths        = NULL) # If the above is filled out
     if (!is.null(TT)) simclass$timegrid <- timegrid(TT, n)
     else simclass$n <- n
 
     class(simclass) <- "simulateclass"

     return(simclass)
}

timegrid <- function(TT, n) {
     s  <- ceiling(n*TT) # total number of time steps
     TT <- s/n           #  = time steps * dt
     
     timegrid <- list(TT = TT,  #  = time steps * dt   maturity (years)
                 n     = n,     # time steps per year
                 dt    = 1/n,   # time step length
                 s     = s,     # total number of time steps
                 t     = seq(0, TT, length = s + 1)) # time step vector
     
     class(timegrid) <- "timegrid"
     return(timegrid)
}

simulategrid <- function(timegrid, t, k) {
     
     simgrid <- volatilitygrid( t = align(t, timegrid$t),
                         k = monotonize(k),
                         impvol = matrix(NA,length(k),length(t)))
     simgrid$s <- index(t,timegrid$t)

     simgrid$siminfo   <- NULL
     
     class(simgrid) <- c("simulategrid","volatilitygrid")
     return(simgrid)
}

subset.simgrid <- function(simgrid, t = simgrid$t, k = simgrid$k) {
     simgrid   <- subset.volatilitygrid(simgrid)
     simgrid$s <- simgrid$s[index(t,simgrid$t)]
     return(simgrid)
}

empiricalgrid <- function(simgrid, volgrid) {

     empgrid <- volatilitygrid( t = align(volgrid$t, simgrid$t),
                         k = align(volgrid$k, simgrid$k),
                         impvol = volgrid$impvol,
                         prices = volgrid$prices)
     
     class(empgrid) <- c("empiricalgrid","volatilitygrid")
     return(empgrid)
}

#########################
##### SET FUNCTIONS #####
#########################

setsimgrid <- function(simclass, t, k) {
     if (is.null(simclass$timegrid)) {
          simclass$timegrid <- timegrid(max(t), simclass$n)
          simclass$n <- NULL
     }
     if (any(t > simclass$timegrid$TT)) warning("some maturities were cut off")
     simclass$simgrid <- simulategrid(simclass$timegrid, t[t <= simclass$timegrid$TT], k)
     return(simclass)
}

setempgrid <- function(simclass, volgrid, calweights = NULL) {
     if (is(volgrid, "dateclass")) volgrid <- volgrid$volgrid
     if (is.null(simclass$simgrid)) simclass <- setsimgrid(simclass, volgrid$t, volgrid$k)
     
     simclass$empgrid <- empiricalgrid(simclass$simgrid, volgrid)
     if (is.null(calweights)) calweights <- matrix(data = 1, nrow = length(volgrid$k), ncol = length(volgrid$t))
     simclass <- set_cal_weights(simclass, calweights)

     return(simclass)
}

changetimegrid <- function(simclass, TT = NULL, n = NULL, reset = FALSE) {
     if (is.null(TT)) TT <- simclass$timegrid$TT
     if (is.null(n)) n   <- simclass$timegrid$n
     simclass$timegrid   <- timegrid(TT,n)
    
     if (reset == FALSE) {
          if (!is.null(simclass$simgrid)) simclass <- setsimgrid(simclass, simclass$simgrid$t, simclass$simgrid$k)
          if (!is.null(simclass$empgrid)) simclass <- setempgrid(simclass, simclass$empgrid)
     } else {
          simclass$simgrid    <- NULL
          simclass$empgrid    <- NULL
          simclass$calweights <- NULL
     }
     simclass$paths <- pathsclass()
     
     return(simclass)
}

set_cal_weights <- function(simclass, weights) {
     if (is.null(simclass$empgrid)) stop("empirical grid must be initialized")
     if (!all(dim(weights) == dim(simclass$empgrid$impvol))) stop("weights must have the right dimensions")

     colnames(weights) <- pround(simclass$empgrid$t)
     rownames(weights) <- pround(simclass$empgrid$k)
     simclass$calweights <- weights
     return(simclass)
}

setweightsto <- function(simclass, values, t = NULL, k = NULL) {
     simclass$calweights[index(k, simclass$empgrid$k), index(t, simclass$empgrid$t)] <- values
     return(simclass)
}

setseed <- function(simclass) {
     if (simclass$seed == -1) {
          simclass$lastsimseed <- -1
     } else {
          if (simclass$seed == 0) simclass$lastsimseed <- floor(as.numeric(Sys.time()))
          else simclass$lastsimseed <- simclass$seed
          set.seed(simclass$lastsimseed)
     }
     
     return(simclass)
}

###################
##### SUMMARY #####
###################

summary.simulateclass <- function(simclass, info = 0) {
     catpas("\ntimegrid variables: TT = ", simclass$timegrid$TT, ", n = ", simclass$timegrid$n, ", N = ", simclass$N, ", seed = ", simclass$seed)
     catpas("\nautomatic variables: dt = ", simclass$timegrid$dt, ", s = ", simclass$timegrid$s)
     
     if (!is.null(simclass$vars)) {
          catpas("\n\nclass variables: ")
          for (i in 1:length(simclass$vars)) {
               catpas(names(simclass$vars)[i], " = ", pround(simclass$vars[[i]]))
               if (i < length(simclass$vars)) catpas(", ")
          }
     }

     catpas("\n\nsimgrid: ", !is.null(simclass$simgrid), " empgrid: ", !is.null(simclass$empgrid))
}

getvars <- function(object, ...) UseMethod("getvars", object)

getvars.simulateclass <- function(simclass, varnames = NULL, digits = NULL) {
     vars <- matrix(as.numeric(simclass$vars), nrow = 1)
     colnames(vars) <- names(simclass$vars)
     if (is.null(varnames)) varnames <- names(simclass$vars)
     if (!is.null(digits)) vars <- round(vars, digits = digits)
     vars[,varnames]
}

########################
##### BENCHMARKING #####
########################

benchmodel <- function(simclass, t = 1, k = 0, n = NULL, times = 10, simfunc = simulate, pricefunc = price, ...) {
     if (is.null(n) && is.null(simclass$timegrid)) n <- simclass$n
     else if(is.null(n)) n <- simclass$timegrid$n
     
     simclass <- changetimegrid(simclass, n = n, TT = max(t), reset = TRUE)
     simclass <- setvars(simclass, ...)
     simclass <- setsimgrid(simclass, t = t, k = k)
     
     benchf <- function() {
          simclass <- simfunc(simclass)
          simclass <- pricefunc(simclass)
     }
     
     mean(microbenchmark(benchf(), times = times)$time)/10^9
}

################
##### PLOT #####
################

plot.simulateclass <- function(simclass, style = c("default","both","simulation","empirical","empirical3d","simulated3d"), mfrow = TRUE, pricetype = c("impvol","prices"), ...) {
     
     if (match.arg(pricetype) == "prices") {
          if (is.null(simclass$simgrid$prices)) simclass$simgrid$prices <- getprices(simclass$simgrid)
          simclass$simgrid$impvol <- simclass$simgrid$prices
          if (is.null(simclass$empgrid$prices)) simclass$empgrid$prices <- getprices(simclass$empgrid)
          simclass$empgrid$impvol <- simclass$empgrid$prices
     }
     
     style <- match.arg(style)
     
     if (style == "default") {
          style <- "both"
          if (is.null(simclass$empgrid)) style <- "simulation"
     }
     
     if (style == "both") {
          if (is.null(simclass$simgrid) || is.null(simclass$empgrid)) stop("Both egrid and ogrid must be initilized")
          
          t <- simclass$empgrid$t
          k <- simclass$empgrid$k
          
          empvol <- simclass$empgrid$impvol
          simvol <- getsmiles(simclass$simgrid, t, k)
          
          if (mfrow) par(mfrow = getmfrow(length(t)))
          
          for (i in 1:length(t)) {
               plot(x = k, y = empvol[,i], ylim = c(min(empvol,simvol),max(empvol,simvol)), xlab = "k", ylab = "impvol", main = paste("T = ", t[i]))
               lines(x = k, y = simvol[,i], col = 2)
          }
          
     } else if (style == "simulation") {
          if (is.null(simclass$simgrid)) stop("simgrid must be initilized")
          plot(simclass$simgrid, mfrow = mfrow, ...)
     } else if (style == "empirical") {
          if (is.null(simclass$empgrid)) stop("empgridmust be initilized")
          plot(simclass$empgrid, mfrow = mfrow,  ...)  
     } else if (style == "empirical3d") {
          if (is.null(simclass$empgrid)) stop("empgridmust be initilized")
          plot_ly(y = simclass$empgrid$t, x = simclass$empgrid$k, z = simclass$empgrid$impvol, type = "surface")
     } else if (style == "simulated3d") {
          if (is.null(simclass$simgrid)) stop("simgrid must be initilized")
          plot_ly(y = simclass$simgrid$t, x = simclass$simgrid$k, z = simclass$simgrid$impvol, type = "surface")
     }
}

plotmultiple <- function(simclasslist) {
     volgridlist <- lapply(simclasslist, function(simclass) simclass$simgrid)
     plotmultiplevolgrids(volgridlist)
}