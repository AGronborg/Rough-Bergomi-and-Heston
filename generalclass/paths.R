###########################
### IMPORTS & LIBRARIES ###
###########################

#############
### CLASS ###
#############
pathsclass <- function() {
     paths <- list()
     class(paths) <- "pathsclass"
     
     return(paths)
}

#########################
##### SET FUNCTIONS #####
#########################

deletepaths <- function(someclass) UseMethod("deletepaths", someclass)

deletepaths.simulateclass <- function(simclass) {
     simclass$paths <- pathsclass()
     return(simclass)
}

################
##### PLOT #####
################

plot.matrix <- function(mat, paths = 1:2, ...) {

     ymin  <- min(mat[paths,])
     ymax  <- max(mat[paths,])
     
     plotfunc <- function(i) {
          if (i == 1) plot(1:ncol(mat), mat[paths[i],], type = "l", col = i, ylim = c(ymin,ymax), ...)   
          else lines(1:ncol(mat), mat[paths[i],], type = "l", col = i)
     }
     
     noreturn <- sapply(1:length(paths), plotfunc)
}

plot.pathsclass <- function(pclass, variables = NULL, paths = 1:2, mfrow = TRUE, ...) {
     if (is.null(variables)) variables <- names(pclass)
     if (mfrow) par(mfrow = getmfrow(length(variables)))
     plotfunc <- function(var) plot(pclass[[var]], paths = paths, main = var)
     noreturn <- sapply(variables, plotfunc)
}