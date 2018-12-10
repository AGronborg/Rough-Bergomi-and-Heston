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

plot.matrix <- function(mat, paths = 1:2, type = "l", ylim = c(min(mat[paths,]),max(mat[paths,])), makeplot = TRUE, col = NULL, ...) {

     plotfunc <- function(i) {
          if (i == 1 && makeplot) plot(1:ncol(mat), mat[paths[i],], type = type, col = ifelse(is.null(col),i,col), ylim = ylim, ...)   
          else lines(1:ncol(mat), mat[paths[i],], type = type, col = ifelse(is.null(col),i,col))
     }
     
     noreturn <- sapply(1:length(paths), plotfunc)
}

plot.pathsclass <- function(pclass, variables = NULL, paths = 1:2, mfrow = TRUE, main = NULL, ...) {
     if (is.null(variables)) variables <- names(pclass)[sapply(pclass, is.matrix)]
     if (mfrow) par(mfrow = getmfrow(length(variables)))
     plotfunc <- function(var) plot(pclass[[var]], paths = paths, main = ifelse(is.null(main),var,main), ...)
     noreturn <- sapply(variables, plotfunc)
}