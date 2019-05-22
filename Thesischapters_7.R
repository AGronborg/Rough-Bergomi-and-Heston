##############################$
##### IMPORTS & LIBRARIES #####
##############################$

### IMPORTS ###
source("thesislib.R")
library(beepr)

#########################################$
##### CHAPTER 7 - VARIANCE REDUCTION #####
#########################################$

##### FUNCTIONS #####

getinfo <- function(classlist, exact = NULL) {
     info <- list()
     
     estimates  <- sapply(classlist, function(x) x$simgrid$impvol)
     times      <- sapply(classlist, simpricetime, units = "sec", digits = 6)
     
     info$k         <- classlist[[1]]$simgrid$k
     info$t         <- classlist[[1]]$simgrid$t
     
     info$estimates <- estimates
     info$means     <- apply(estimates, 1, mean)
     info$vars      <- apply(estimates, 1, var)
     info$stddev    <- apply(estimates, 1, sd)
     info$times     <- times
     info$tau       <- mean(times)
     
     if (is.null(exact)) exact <- info$means
     if (is(exact, "volatilitygrid")) exact <- exact$impvol
     
     info$exact  <- exact
     info$bias   <- as.numeric(info$means - exact)
     
     N           <- length(classlist)
     strikes     <- length(classlist[[1]]$simgrid$k)
     bvar        <- function(estimates, exact, N) 1/(N-1) * sum((estimates - exact)^2)
     info$bvars  <- sapply(1:strikes, function(i) bvar(estimates[i,], exact[i], N))
     
     
     info$phi    <- sqrt(         mean(info$bvars))
     info$theta  <- sqrt(info$tau*mean(info$bvars))
     info$phi2   <-               mean(info$bvars)
     info$theta2 <-      info$tau*mean(info$bvars)
     
     return(info)
}

plotlist <- function(classlist, exact = NULL, xlim = c(-0.02,0.02), breaks = seq(-1, 1, length.out = 1000), ...) {
     info      <- getinfo(classlist, exact)
     estimates <- info$estimates
     
     for (i in 1:length(info$k)) {
          k <- info$k[i]
          
          hist(estimates[i,]-info$exact[i], main = paste("k = ", round(k,2)), probability = TRUE, xlim = xlim, breaks = breaks, ...)
          
          x <- seq(xlim[1], xlim[2], length.out = 100)
          y <- dnorm(x = x, mean = info$bias[i], sd = info$stddev[i])
          lines(x, y, col = 2, lwd = 2)
          abline(v = info$bias[i], col = 4, lwd = 2, lty = 2)
     }
}

gettablelist <- function(classlistlist, exact = NULL, names = "") {
     
     getrow <- function(classlist) {
          info <- getinfo(classlist, exact)
          
          row <- c(info$stddev, info$tau, info$phi*100, info$theta*100)
          
          return(row)
     }
     
     k     <- getinfo(classlistlist[[1]])$k
     table <- t( sapply(classlistlist, getrow) )
     table <- cbind(table, table[1,ncol(table)]/table[,ncol(table)])
     
     colnames(table) <- c(round(k,2), "tau", "phi", "theta", "imp")
     rownames(table) <- names
     return(table)
}

#######################################
##### CONVERGENCE (N) - 3 strikes #####
#######################################

# Strikes
strikes <- c(-0.1787,0.0000,0.1041)

# Exact
exact1  <- c(0.2961,0.2061,0.1576)

# Base
rbclasses11 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312*4, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = 0.25, strikes)
     rbclass <- simulate_rb_standard(rbclass)
     rbclass <- price_standard(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses11[[N]] <- rbclass
}

save(rbclasses11, file = "data/cht7_volgrids_base.RData")

# Antithetic (Base)
rbclasses12 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312*4, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = 0.25, strikes)
     rbclass <- simulate_rb(rbclass, skip = "S1", antithetic = TRUE)
     rbclass <- price_standard(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses12[[N]] <- rbclass
}

save(rbclasses12, file = "data/cht7_volgrids_antithetic.RData")

# Conditional
rbclasses13 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312*4, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = 0.25, strikes)
     rbclass <- simulate_rb(rbclass)
     rbclass <- price_rb_mixed(rbclass, estimator = "conditional")
     rbclass <- deletepaths(rbclass)
     rbclasses13[[N]] <- rbclass
}

save(rbclasses13, file = "data/cht7_volgrids_conditional.RData")

# Controlled
rbclasses14 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312*4, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = 0.25, strikes)
     rbclass <- simulate_rb_standard(rbclass)
     rbclass <- price_rb_mixed(rbclass, estimator = "controlled")
     rbclass <- deletepaths(rbclass)
     rbclasses14[[N]] <- rbclass
}

save(rbclasses14, file = "data/cht7_volgrids_controlled.RData")

# Mixed (no antithetic)
rbclasses15 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312*4, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = 0.25, strikes)
     rbclass <- simulate_rb_mixed(rbclass)
     rbclass <- price_rb_mixed(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses15[[N]] <- rbclass
}

save(rbclasses15, file = "data/cht7_volgrids_mixed.RData")

# Mixed (antithetic)
rbclasses16 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312*4, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = 0.25, strikes)
     rbclass <- simulate_rb_antimixed(rbclass)
     rbclass <- price_rb_mixed(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses16[[N]] <- rbclass
}

save(rbclasses16, file = "data/cht7_volgrids_antimixed.RData")

# Load
load(file = "data/cht7_volgrids_base.RData")
load(file = "data/cht7_volgrids_antithetic.RData")
load(file = "data/cht7_volgrids_conditional.RData")
load(file = "data/cht7_volgrids_controlled.RData")
load(file = "data/cht7_volgrids_mixed.RData")
load(file = "data/cht7_volgrids_antimixed.RData")

rbclasseslist <- list(rbclasses11,rbclasses12,rbclasses13,rbclasses14,rbclasses15,rbclasses16)
table <- gettablelist(rbclasseslist, exact = c(0.2961,0.2061,0.1576), names = c("base", "antithetic", "conditional", "controlled", "mixed", "antimixed"))
round(table, 3)

par(mfrow = c(6,3), mar = c(2,2,1,0.5))
lapply(rbclasseslist, plotlist, exact = c(0.2961,0.2061,0.1576), xlim = c(-0.06,0.06))

# Table
xtable(table, digits = 3)

# Plot
pdf("../Thesis LaTeX/Images/cht7_plot3.pdf")
par(mfrow = c(6,3), mar = c(2,2,1,0.5))
lapply(rbclasseslist, plotlist, exact = c(0.2961,0.2061,0.1576), xlim = c(-0.06,0.06))
dev.off()

#########################################
##### CONVERGENCE (N) - FULL SMILES #####
#########################################

# Strikes 
strikes <- loaddate(date = "2006-10-30")$volgrid$k
t       <- 1

# Exact
rbclass <- roughbergomiclass(n = 312/t, N = 20000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
rbclass <- setsimgrid(rbclass, t = t, k = strikes)
exact2  <- splitprice(rbclass, times = 20, sim_func = function(x) simulate_rb_standard(x, kappa = 1, b = "optimal"), price_func = function(...) price_standard(...))

save(exact2, file = "data/cht7_volgrids3_exact.RData")
load(file = "data/cht7_volgrids3_exact.RData")

# Base
rbclasses21 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312/t, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = t, strikes)
     rbclass <- simulate_rb_standard(rbclass)
     rbclass <- price_standard(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses21[[N]] <- rbclass
}

save(rbclasses21, file = "data/cht7_volgrids3_base.RData")

# Antithetic (Base)
rbclasses22 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312/t, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = t, strikes)
     rbclass <- simulate_rb(rbclass, skip = "S1", antithetic = TRUE)
     rbclass <- price_standard(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses22[[N]] <- rbclass
}

save(rbclasses22, file = "data/cht7_volgrids3_antithetic.RData")

# Conditional
rbclasses23 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312/t, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = t, strikes)
     rbclass <- simulate_rb(rbclass)
     rbclass <- price_rb_mixed(rbclass, estimator = "conditional")
     rbclass <- deletepaths(rbclass)
     rbclasses23[[N]] <- rbclass
}

save(rbclasses23, file = "data/cht7_volgrids3_conditional.RData")

# Controlled
rbclasses24 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312/t, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = t, strikes)
     rbclass <- simulate_rb_standard(rbclass)
     rbclass <- price_rb_mixed(rbclass, estimator = "controlled")
     rbclass <- deletepaths(rbclass)
     rbclasses24[[N]] <- rbclass
}

save(rbclasses24, file = "data/cht7_volgrids3_controlled.RData")

# Mixed (no antithetic)
rbclasses25 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312/t, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = t, strikes)
     rbclass <- simulate_rb_mixed(rbclass)
     rbclass <- price_rb_mixed(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses25[[N]] <- rbclass
}

save(rbclasses25, file = "data/cht7_volgrids3_mixed.RData")

# Mixed (antithetic)
rbclasses26 <- list()

for (N in 1:1000) {
     print(N)
     rbclass <- roughbergomiclass(n = 312/t, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
     rbclass <- setsimgrid(rbclass, t = t, strikes)
     rbclass <- simulate_rb_antimixed(rbclass)
     rbclass <- price_rb_mixed(rbclass)
     rbclass <- deletepaths(rbclass)
     rbclasses26[[N]] <- rbclass
}

save(rbclasses26, file = "data/cht7_volgrids3_antimixed.RData")

# Load (t=0.25)
# load(file = "data/cht7_volgrids2_base.RData")
# load(file = "data/cht7_volgrids2_antithetic.RData")
# load(file = "data/cht7_volgrids2_conditional.RData")
# load(file = "data/cht7_volgrids2_controlled.RData")
# load(file = "data/cht7_volgrids2_mixed.RData")
load(file = "data/cht7_volgrids2_antimixed.RData")
load(file = "data/cht7_volgrids2_exact.RData")

# Load (t=1)
load(file = "data/cht7_volgrids3_base.RData")
load(file = "data/cht7_volgrids3_antithetic.RData")
load(file = "data/cht7_volgrids3_conditional.RData")
load(file = "data/cht7_volgrids3_controlled.RData")
load(file = "data/cht7_volgrids3_mixed.RData")
load(file = "data/cht7_volgrids3_antimixed.RData")
load(file = "data/cht7_volgrids3_exact.RData")

rbclasseslist <- list(rbclasses21,rbclasses22,rbclasses23,rbclasses24,rbclasses25,rbclasses26)
table <- gettablelist(rbclasseslist, exact = exact2, names = c("base", "antithetic", "conditional", "controlled", "mixed", "antimixed"))

round(table, 3)

par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses21, exact = exact2, xlim = c(-0.2,0.2))
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses22, exact = exact2, xlim = c(-0.2,0.2))
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses23, exact = exact2, xlim = c(-0.2,0.2))
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses24, exact = exact2, xlim = c(-0.05,0.05))
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses25, exact = exact2, xlim = c(-0.05,0.05))
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses26, exact = exact2, xlim = c(-0.02,0.02))

# table
xtable(table, digits = 3)

# plots
pdf("../Thesis LaTeX/Images/cht7_plot4.pdf")
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses21, exact = exact2, xlim = c(-0.2,0.2))
dev.off()

pdf("../Thesis LaTeX/Images/cht7_plot5.pdf")
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses26, exact = exact2, xlim = c(-0.02,0.02))
dev.off()
