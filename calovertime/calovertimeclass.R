###########################
### IMPORTS & LIBRARIES ###
###########################

source("generalclass/simclass.R")
source("rbclass/rbclass.R")
source("hestonclass/hestonclass.R")

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
options(xtable.comment = FALSE)

###############
### CLASSES ###
###############
calovertimeclass <- function(datesclass, simclass) {
     cotclass <- list(dates        = onlydates(datesclass),
                      mainsimclass = simclass,
                      simclasses   = lapply(1:datesclass$info$nrobs, function(x) return(simclass)),
                      info         = datesclass$info)
     class(cotclass) <- "calovertimeclass"
     
     cotclass <- setempgrid_all(cotclass, datesclass)

     return(cotclass)
}

#########################
##### SET FUNCTIONS #####
#########################

runonall <- function(cotclass, func, ...) {
     for (i in 1:cotclass$info$nrobs) {
          cotclass$simclasses[[i]] <- func(cotclass$simclasses[[i]], ...)
     }
     return(cotclass)
}

setsimgrid_all <- function(cotclass, datesclass) {
     for (i in 1:cotclass$info$nrobs) {
          cotclass$simclasses[[i]] <- setsimgrid(cotclass$simclasses[[i]], datesclass[[i]]$volgrid$t, datesclass[[i]]$volgrid$k)
     }
     return(cotclass)
}

setempgrid_all <- function(cotclass, datesclass) {
     for (i in 1:cotclass$info$nrobs) cotclass$simclasses[[i]] <- setempgrid(cotclass$simclasses[[i]], datesclass[[i]])
     return(cotclass)
}

setvars.calovertimeclass <- function(cotclass, ...) runonall(cotclass, setvars, ...)

setvars_tofirst <- function(cotclass) {
     for (i in 2:cotclass$info$nrobs) 
          cotclass$simclasses[[i]]$vars <- cotclass$simclasses[[1]]$vars
     return(cotclass)
}

deletepaths.calovertimeclass <- function(cotclass) runonall(cotclass, deletepaths)

#########################
##### GET FUNCTIONS #####
#########################

getvars.calovertimeclass <- function(cotclass, varnames = NULL, dates = NULL, digits = NULL) {
     
     nrobs  <- cotclass$info$nrobs
     nrvars <- length(cotclass$simclasses[[1]]$vars)
     
     vars <- matrix(NA, ncol = nrvars, nrow = nrobs)
     rownames(vars) <- showdates(cotclass$dates)
     colnames(vars) <- names(cotclass$simclasses[[1]]$vars)
     
     for (i in 1:nrobs) vars[i, ] <- as.numeric(cotclass$simclasses[[i]]$vars)
     
     if (is.null(varnames)) varnames <- names(cotclass$simclasses[[1]]$vars)
     if (is.null(dates))    dates    <- 1:nrobs
     if (is.numeric(dates)) dates    <- showdates(cotclass$dates)[dates]
     
     if (!is.null(digits)) vars <- round(vars, digits = digits)
     
     vars[dates,varnames]
}

#####################
##### CALIBRATE #####
#####################

calibrate_first <- function(cotclass, ...) {
     cotclass$info$firstcalstarttime <- Sys.time()
     cotclass$simclasses[[1]] <- calibrate(cotclass$simclasses[[1]], ...)
     cotclass$info$firstcalendtime <- Sys.time()
     return(cotclass)
}

calibrate.calovertimeclass <- function(cotclass, skipfirst = FALSE, delpaths = TRUE, ...) {
     cotclass$info$calstarttime <- Sys.time()
     for (i in (1 + skipfirst):cotclass$info$nrobs) {
          cotclass$simclasses[[i]] <- calibrate(cotclass$simclasses[[i]], ...)
          if (delpaths) cotclass$simclasses[[i]] <- deletepaths(cotclass$simclasses[[i]])
     }
     cotclass$info$calendtime <- Sys.time()
     return(cotclass)
}

caltime.calovertimeclass <- function(cotclass, units = "secs", digits = 0) {
     round(sum( sapply(cotclass$simclasses, caltime, units = units) ), digits = digits)
     #caltime <- as.numeric(difftime(cotclass$info$calendtime,cotclass$info$calstarttime, units = units))
     #firstcaltime <- as.numeric(difftime(cotclass$info$firstcalendtime,cotclass$info$firstcalstarttime, units = units))
     #return(round(firstcaltime+caltime, digits = digits))
}

caltimes <- function(cotclass, units = "secs", digits = 0) {
     caltimes <- matrix(sapply(cotclass$simclasses, caltime, units = units), ncol = 1)
     colnames(caltimes) <- units
     rownames(caltimes) <- showdates(cotclass$dates)
     return(caltimes)
}

##################
##### XTABLE #####
##################

xtable_h_rb <- function(hcotclass, rbcotclass, pars = rep(2, 9), dates = 1:hcotclass$info$nrobs, ...) {
     hvars <- getvars(hcotclass)[dates,which(pars[1:5] > 0)]
     rbvars <- getvars(rbcotclass)[dates,which(pars[6:9] > 0)]
     
     hvars[2:length(dates),which(pars[1:5] == 1)]  <- NA
     rbvars[2:length(dates),which(pars[6:9] == 1)] <- NA
     
     vars <- cbind(hvars, rep(NA, length(dates)), rbvars)

     xt <- xtable(vars, ...)
     names(xt) <- c("$\\lambda$", "$\\bar{\\nu}$", "$\\nu_0$", "$\\eta$", "$\\rho$", "","$\\alpha$", "$\\rho$", "$\\eta$", "$\\xi$")
     align(xt) <- rep("c", 11)
     
     print(xt, sanitize.colnames.function = identity)
}

###################
##### SUMMARY #####
###################

summary.calovertimeclass <- function(cotclass) {
     vars     <- pround(getvars(cotclass))
     
     errors   <- pround(sapply(cotclass$simclasses, function(x) x$calinfo$op$value))
     vars     <- cbind(vars, errors)
     colnames(vars)[ncol(vars)] <- "sse"
     
     ctimes   <- caltimes(cotclass)
     vars     <- cbind(vars, ctimes)
     
     return(vars)
}

################
##### PLOT #####
################

plot.calovertimeclass <- function (cotclass, ylim = c(0,0.8)) multiplot_calovertime(list(cotclass), ylim = ylim)

multiplot_calovertime <- function(cotclasslist, t = NULL, ...) {
     if (is.null(t) || length(t) > 1) multiplot_calovertime_mt(cotclasslist, t, ...)
     else if (length(t) == 1) multiplot_calovertime_1t(cotclasslist, t, ...)
}

multiplot_calovertime_1t <- function(cotclasslist, t, withempirical = TRUE, mfrow = TRUE, ylim = NULL) {

     s <- index(t, cotclasslist[[1]]$simclasses[[1]]$simgrid$t)

     if (is.null(ylim)) {
          ymin <- apply(sapply(cotclasslist, function(x) sapply(x$simclasses, function(y) min(y$simgrid$impvol[,s])) ), 1, min)
          ymax <- apply(sapply(cotclasslist, function(x) sapply(x$simclasses, function(y) max(y$simgrid$impvol[,s])) ), 1, max)
     } else {
          ymin <- rep(ylim[1], cotclasslist[[1]]$info$nrobs)
          ymax <- rep(ylim[2], cotclasslist[[1]]$info$nrobs)
     }
     
     if (mfrow) par(mfrow = getmfrow(cotclasslist[[1]]$info$nrobs))
     for (i in 1:cotclasslist[[1]]$info$nrobs) {
          first <- TRUE
          gridnum   <- 1
          for (cotclass in cotclasslist) {
               volgrid <- cotclass$simclasses[[i]]$simgrid
               
               if (first && withempirical) {
                    y <- cotclass$simclasses[[i]]$empgrid$impvol[,s]
                    plot(x = cotclass$simclasses[[1]]$empgrid$k, y = y, 
                               type = "p", xlab = "k", ylab = "impvol", col = gridnum,
                               ylim = c(min(ymin[i],y),max(ymax[i],y)), #xlim = c(xmin,xmax), 
                               main = paste(cotclass$dates[[i]]$date, " (T = ", round(volgrid$t[s],2), ")", sep = ""))
               }
               
               if (first && !withempirical) plot(x = volgrid$k, y = volgrid$impvol[,s], 
                                   type = "l", xlab = "k", ylab = "impvol", col = gridnum,
                                   ylim = c(ymin[i],ymax[i]), #xlim = c(xmin,xmax), 
                                   main = paste(cotclass$dates[[i]]$date, " (T = ", round(volgrid$t[s],2), ")", sep = ""))
               else lines(x = volgrid$k, y = volgrid$impvol[,s], col = gridnum)
               gridnum <- gridnum + 1
               first   <- FALSE
          }
     }
}

multiplot_calovertime_mt <- function(cotclasslist, t = NULL, withempirical = TRUE, mfrow = TRUE, ylim = NULL) {
     
     if (is.null(t)) t <- cotclasslist[[1]]$simclasses[[1]]$simgrid$t
     
     s <- rep(NA, length(t))
     for (i in 1:length(t)) s[i] <- index(t[i], cotclasslist[[1]]$simclasses[[1]]$simgrid$t)
     
     if (is.null(ylim)) {
          #ymin <- apply(sapply(cotclasslist, function(x) sapply(x$simclasses, function(y) min(y$simgrid$impvol[,s])) ), 1, min)
          #ymax <- apply(sapply(cotclasslist, function(x) sapply(x$simclasses, function(y) max(y$simgrid$impvol[,s])) ), 1, max)
     } else {
          ymin <- rep(ylim[1], cotclasslist[[1]]$info$nrobs)
          ymax <- rep(ylim[2], cotclasslist[[1]]$info$nrobs)
     }
     
     #if (mfrow) par(mfrow = getmfrow(cotclasslist[[1]]$info$nrobs*length(t)))
     if (mfrow) par(mfrow = c(length(t), cotclasslist[[1]]$info$nrobs))
     for (j in 1:length(s)) {
          for (i in 1:cotclasslist[[1]]$info$nrobs) {
               first <- TRUE
               gridnum   <- 1
               for (cotclass in cotclasslist) {
                    volgrid <- cotclass$simclasses[[i]]$simgrid
                    
                    #if (j == 1) main <- paste(cotclass$dates[[i]]$date, " (T = ", volgrid$t[s[j]], ")")
                    #else main <- "" 
                    
                    main <- paste(cotclass$dates[[i]]$date, " (T = ", round(volgrid$t[s[j]],2), ")")
                    
                    if (first && withempirical) {
                         y <- cotclass$simclasses[[i]]$empgrid$impvol[,s[j]]
                         plot(x = cotclass$simclasses[[1]]$empgrid$k, y = y, 
                              type = "p", xlab = "k", ylab = "impvol", col = gridnum,
                              ylim = c(min(ymin[i],y),max(ymax[i],y)), #xlim = c(xmin,xmax), 
                              main = main)
                    }
                    
                    if (first && !withempirical) plot(x = volgrid$k, y = volgrid$impvol[,s[j]], 
                                                      type = "l", xlab = "k", ylab = "impvol", col = gridnum,
                                                      ylim = c(ymin[i],ymax[i]), #xlim = c(xmin,xmax), 
                                                      main = main)
                    else lines(x = volgrid$k, y = volgrid$impvol[,s[j]], col = gridnum + 1)
                    gridnum <- gridnum + 1
                    first   <- FALSE
               }
          }
     }
}
