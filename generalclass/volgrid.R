###########################
### IMPORTS & LIBRARIES ###
###########################

###############
### CLASSES ###
###############

volatilitygrid <- function(t, k, impvol, prices = NULL) {
     
     if (!is.matrix(impvol)) stop("implied volatility must be a matrix")
     rownames(impvol) <- NULL
     colnames(impvol) <- NULL
     
     if (!(ncol(impvol) == length(t)) || !(nrow(impvol) == length(k)) )
          stop("t, k and impvol must have the right dimensions (t = ncol(impvol))")
     
     volgrid <- list(t = t, k = k, impvol = impvol, prices = prices) # k = log(K)
     class(volgrid) <- "volatilitygrid"
     
     return(volgrid)
}

############################
##### CALCULATE PRICES #####
############################

getprices <- function(volgrid) {
     t      <- volgrid$t
     k      <- volgrid$k
     impvol <- volgrid$impvol
     
     prices <- matrix(NA, nrow = nrow(impvol), ncol = ncol(impvol))
     for (i in 1:length(t)) prices[,i] <- vec_bs(Fwd = 1, K = exp(k), V = impvol[,i]^2*t[i])
     
     return(prices)
}

getimpvol <- function(volgrid) {
     t      <- volgrid$t
     k      <- volgrid$k
     prices <- volgrid$prices
     
     impvol <- matrix(NA, nrow = nrow(prices), ncol = ncol(prices))
     for (i in 1:length(t)) impvol[,i] <- vec_bsinv(P = prices[,i], Fwd = 1, K = exp(k), TT = t[i])
     
     return(impvol)
}

#############################
##### VOLGRID FUNCTIONS #####
#############################

subset.volatilitygrid <- function(volgrid, t = NULL, k = NULL) {
     
     if (is.null(t)) t <- volgrid$t
     if (is.null(k)) k <- volgrid$k
     
     t <- align(t,volgrid$t)
     k <- align(k,volgrid$k)

     volgrid$t <- t
     volgrid$k <- k
     volgrid$impvol <- as.matrix(volgrid$impvol[index(k,volgrid$k),index(t,volgrid$t)])
     if (!is.null(volgrid$prices)) volgrid$prices <- as.matrix(volgrid$prices[index(k,volgrid$k),index(t,volgrid$t)])
     
     return(volgrid)
}

getsmiles <- function(volgrid, t = NULL, k = NULL, pricetype = c("impvol","prices")) {
     if (match.arg(pricetype) == "impvol") return(subset(volgrid, t, k)$impvol)
     else if (match.arg(pricetype) == "prices") return(subset(volgrid, t, k)$prices)
}

align_grid <- function(volgrid1,volgrid2) {
     volgrid1$t <- align(volgrid1$t,volgrid2$t)
     volgrid1$k <- align(volgrid1$k,volgrid2$k)
     volatilitygrid(volgrid1$t,volgrid1$k,volgrid1$impvol)
}

stdvolgrid <- function(t, k) {
     if (is.number(t)) t <- seq(0,1,length.out = t)
     if (is.number(k)) k <- seq(-0.3,0.3,length.out = k)
     volatilitygrid(t, k, matrix(sample(1:40*0.01, length(t)*length(k), TRUE), nrow = length(k), ncol = length(t)))
}


################
##### PLOT #####
################

plot.volatilitygrid <- function(volgrid, t = NULL, plot3d = FALSE, ...) {
     if (plot3d == TRUE) plot_ly(y = volgrid$t, x = volgrid$k, z = volgrid$impvol, type = "surface")
     else volgrid_multiplot(volgrid, t, ...)
}

volgrid_multiplot <- function(volgrid, t = NULL, mfrow = TRUE, pricetype = c("impvol","prices"), ...) {
     
     if (match.arg(pricetype) == "prices") {
               if (is.null(volgrid$prices)) volgrid$prices <- getprices(volgrid)
               volgrid$impvol <- volgrid$prices
     }
     
     s <- NULL
     if (!is.null(t)) s <- index(t,volgrid$t)
     else s <- 1:length(volgrid$t)
     
     if (mfrow) par(mfrow = getmfrow(length(s)))
     for (i in s) plot(x = volgrid$k, y = volgrid$impvol[,i], xlab = "k", ylab = "impvol", main = paste("T = ", volgrid$t[i]), ...)
}

plotmultiplevolgrids <- function(volgridlist, mfrow = TRUE, pricetype = c("impvol","prices")) {
     
     if (match.arg(pricetype) == "prices") {
          for (i in 1:length(volgridlist)) {
               if (is.null(volgridlist[[i]]$prices)) volgridlist[[i]]$prices <- getprices(volgridlist[[i]])
               volgridlist[[i]]$impvol <- volgridlist[[i]]$prices
          }
     }
     
     xmin <- min(sapply(volgridlist, function(volgrid) min(volgrid$k)))
     xmax <- max(sapply(volgridlist, function(volgrid) max(volgrid$k)))
     
     ymin <- apply(as.matrix(sapply(volgridlist, function(volgrid) apply(volgrid$impvol, 2, min))), 1, min)
     ymax <- apply(as.matrix(sapply(volgridlist, function(volgrid) apply(volgrid$impvol, 2, max))), 1, max)
     
     t <- volgridlist[[1]]$t
     s <- 1:length(t)
     
     if (mfrow) par(mfrow = getmfrow(length(s)))
     for (i in s) {
          first <- TRUE
          gridnum   <- 1
          for (volgrid in volgridlist) {
               if (first) plot(x = volgrid$k, y = volgrid$impvol[,i], 
                               type = "l", xlab = "k", ylab = "impvol",
                               xlim = c(xmin,xmax), ylim = c(ymin[i],ymax[i]),
                               main = paste("T = ", volgrid$t[i]))
               else lines(x = volgrid$k, y = volgrid$impvol[,i], col = gridnum)
               gridnum <- gridnum + 1
               first   <- FALSE
          }
     }
}
