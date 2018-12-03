####################
##### SIMULATE #####
####################

simulate               <- function(simclass, ...) UseMethod("simulate", simclass)
simulate.simulateclass <- function(simclass) return(simclass)

##############################
##### SIMULATE FUNCTIONS #####
##############################

simtime <- function(object, ...) UseMethod("simtime", object)
simtime.simulateclass <- function(simclass, units = "auto", digits = 0) gettime(simclass$siminfo, units, digits)