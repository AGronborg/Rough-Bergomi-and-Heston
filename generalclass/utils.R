###########################
### IMPORTS & LIBRARIES ###
###########################

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
options(xtable.comment = FALSE)

#########################
### UTILITY FUNCTIONS ###
#########################

pround <- function(x) round(x, 3)

catpas <- function(...) {
     txt <- paste(..., sep = "")
     cat(txt)
}

is.number <- function(x) {
     is.numeric(x) && length(x) == 1 && !is.nan(x)
}

monotonize <- function(t) {
     if (!all(sort(t) == t)) {
          warning("variable has been sorted")
          t <- sort(t)
     }
     
     if (!(length(unique(t)) == length(t))) {
          warning(paste("variable has been made unique. Length has gone from", length(t), "to", length(unique(t))))
          t <- unique(t)
     }
     
     return(t)
}

align <- function(t1, t2) {
     # if (!all(is.in(t1,t2))) warning("t1 has been rounded off")
     monotonize( t2[which.closest(t1,t2)] )
}

which.closest <- function(t1,t2) {
     sapply(t1, function(x) which.min(abs(x-t2))  )
}

is.in <- function(t1,t2) {
     sapply(t1, function(x) sum(abs(x - t2) < 1e-8) >= 1 ) # tol = 1e-8
}

samesize <- function(x, value = NA) {
     if (is.matrix(x)) {
          matrix(value, nrow = nrow(x), ncol = ncol(x))
     } else if (is.numeric(x)) {
          rep(value, length(x))
     }
}

stdmatrix <- function(rows, cols, ...) {
     matrix(sample(1:40*0.01, rows*cols, TRUE), nrow = rows, ncol = cols, ...)
}

index <- function(t1,t2) {
     t1 <- align(t1,t2)
     sapply( t1, function(x) match(x,t2) )
}

getmfrow <- function(n) {
     if (n > 16) n <- 16
     switch(n,
            c(1,1),c(2,1),c(2,2),c(2,2),
            c(3,2),c(3,2),c(4,2),c(4,2),
            c(3,3),c(4,3),c(4,3),c(4,3),
            c(4,4),c(4,4),c(4,4),c(4,4))
}

se  <- function(x) sqrt(var(x)/length(x))
rse <- function(x) abs(se(x)/mean(x))*100

##### GETTIME #####
gettime <- function(timelist, units = "auto", digits) {
     if (units == "auto") difftime(timelist$endtime,timelist$starttime, units = units)
     else round(as.numeric(difftime(timelist$endtime,timelist$starttime, units = units)), digits = digits)
}

##### Black-Scholes #####
bs <- function(Fwd, K, V, o = "call") {
     w <- 1
     if (o == "put") w <- -1
     else if (o == "otm") w <- 2 * (K > 1) - 1
     else if (o == "itm") w <- 2 * (K < 1) - 1
     
     sv <- sqrt(V)
     d1 <- log(Fwd/K) / sv + 0.5 * sv
     d2 <- d1 - sv
     P  <- w * Fwd * pnorm(w * d1) - w * K * pnorm(w * d2)
     
     return (P)
}
vec_bs <- Vectorize(bs)

##### Black-Scholes implied volatility #####
bsinv <- function(P, Fwd, K, TT, o = "call") {
     w <- 1
     if (o == "put") w <- -1
     else if (o == "otm") w <- 2 * (K > 1) - 1 # Choose call/put whatever is OTM
     else if (o == "itm") w <- 2 * (K < 1) - 1
     
     P <- max(P, max(w * (Fwd - K),0) )
     
     if ( (bs(Fwd, K, 1e-9^2*TT, o) - P)*(bs(Fwd, K, 1e+9^2*TT, o) - P) > 0 ) return(0) # If no solution found
     else impV <- uniroot(function(s) { bs(Fwd, K, s^2*TT, o) - P }, c(1e-9,1e+9))
     
     return (impV$root)
}
vec_bsinv <- Vectorize(bsinv)

##### strikes and deltas #####
delta_put <- function(k,v,S0 = 1) {
     d1 <- (log(S0)-k) / v + 0.5 * v
     pnorm(-d1)
}

delta_call <- function(k,v,S0 = 1) {
     d1 <- (log(S0)-k) / v + 0.5 * v
     pnorm(d1)
}


##################
##### XTABLE #####
##################

xtabletofile <- function(table, filename = "table", path = "tables/") {
     if (filename == "table") {
          i <- 1
          while (file.exists(paste(path, "table", i, ".txt", sep = ""))) i <- i + 1
          filename <- paste("table", i, ".txt", sep = "")
     }
     print(table, file = paste(path, filename, sep = ""))
}
