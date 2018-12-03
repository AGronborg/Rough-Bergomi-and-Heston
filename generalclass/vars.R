###########################
### IMPORTS & LIBRARIES ###
###########################

#############
### CLASS ###
#############
varclass <- function(names, vars) {
     vclass <- list()
     for (i in 1:length(names)) vclass[[names[i]]] <- vars[i]
     class(vclass) <- "varclass"
     
     return(vclass)
}

#########################
##### SET FUNCTIONS #####
#########################

setvars <- function(object, ...) UseMethod("setvars", object)

setvars.varclass <- function(vclass, names = NULL, values = NULL, ...) {
     if (length(names) > 0) 
          for (i in 1:length(names)) 
               vclass[[names[i]]] <- values[i]
          
     varlist <- list(...)
     if (length(varlist) > 0) 
          for (i in 1:length(varlist)) 
               vclass[[names(varlist)[i]]] <- varlist[[i]]
     
     return(vclass)
}

setvars.simulateclass <- function(simclass, N = NULL, seed = NULL, n = NULL, TT = NULL, ...) {
     if(!is.null(N)) {
          simclass$N <- N
          simclass$paths <- pathsclass()
     }
     if(!is.null(seed)) simclass$seed <- seed
     if (!is.null(n) || !is.null(TT)) stop("Use changetimegrid function instead.")
     
     simclass$vars <- setvars(simclass$vars, ...)
     
     return(simclass)
}

as.varclass <- function(object) UseMethod("as.varclass", object)
as.varclass.list <- function(somelist) varclass(names(somelist), as.numeric(somelist))

identifyvars <- function(varnames, ...) {
     possiblevars <- list(...)
     vars <- list()
     
     if (length(varnames) > 0) {
          for (i in 1:length(varnames)) {
               index <- match(varnames[i],names(possiblevars))
               if (!is.na(index)) vars[[varnames[i]]] <- possiblevars[[index]]
          }
     }
     
     return(as.varclass(vars))
}

isvar <- function(vclass, varnames) varnames %in% names(varclass)

getbounds <- function(simclass, varnames) {
     bounds <- list(lb = c(), ub = c())
     for (i in 1:length(varnames)) {
          index <- match(varnames[i],names(simclass$vars))
          if (!is.na(index)) {
               bounds$lb <- c(bounds$lb,simclass$varbounds$lb[index])
               bounds$ub <- c(bounds$ub,simclass$varbounds$ub[index])
          }
     }
     return(bounds)
}

