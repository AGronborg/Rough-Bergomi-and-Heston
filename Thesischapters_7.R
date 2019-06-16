##############################$
##### IMPORTS & LIBRARIES #####
##############################$

### IMPORTS ###
source("thesislib.R")
library(beepr)

#########################################$
##### CHAPTER 7 - VARIANCE REDUCTION #####
#########################################$

### ANTITHETIC PATHS###
rbclass   <- roughbergomiclass(n = 312, N = 2000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = 123)
rbclass   <- setsimgrid(rbclass, t = 1, k = 0)
rbclass   <- simulate_rb(rbclass, antithetic = TRUE)
rbclass   <- rb_add_antithetic_paths(rbclass)

cor(rbclass$paths$S[1:1000,313],rbclass$paths$S[1001:2000,313])

##### CONVERGENCE (n) #####

# For 3 strikes
rb <- roughbergomiclass(N = 200000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = 123)
simulate.roughbergomiclass  <- function(...) simulate_rb_antimixed(...)
price.roughbergomiclass     <- function(...) price_rb_mixed(...)

strikes <- c(-0.1787, 0, 0.1041)
x <- c(1:10 * 40, 6:16 * 80)
y <- sapply(x, function (n) {print(n); fastprice(rb, t = 0.25, k = strikes, n = n, pricetype = "impvol", seed = 123)} )

# save(strikes, x, y, file = "data/cht7_converge_n.RData")
# load(file = "data/cht7_converge_n.RData")

# Absolute error
conf <- 1/sqrt(200) * c(0.0055,0.0027,0.0026) * qnorm(0.975) # 95% confidence interval (error)

pdf("Images/cht7_plot1.pdf")
par(mfrow=c(1,1), mar = c(4.5,4.5,1,0.5))
plot( x/4, abs(y[1,]-0.2961), type = "l", col = 1, ylim = c(0,0.02), xlab = "time steps", ylab = "absolute error")
lines(x/4, abs(y[2,]-0.2061), type = "l", col = 2)
lines(x/4, abs(y[3,]-0.1576), type = "l", col = 3)
abline(h = conf[1:3], col = 1:3, lwd = 2)
dev.off()

# Relative error
conf <- 1/sqrt(200) * c(0.0055,0.0027,0.0026) * qnorm(0.975) / c(0.2961,0.2061,0.1576) * 100 # 95% confidence interval (error)

pdf("Images/cht7_plot2.pdf")
par(mfrow=c(1,1), mar = c(4.5,4.5,1,0.5))
plot( x/4, abs(y[1,]-0.2961)/0.2961*100, type = "l", lwd = 1, col = 1, ylim = c(0,4), xlab = "time steps (n/4)", ylab = "relative absolute error")
lines(x/4, abs(y[2,]-0.2061)/0.2061*100, type = "l", lwd = 1, col = 2)
lines(x/4, abs(y[3,]-0.1576)/0.1576*100, type = "l", lwd = 1, col = 3)
abline(h = conf[1:3], col = 1:3, lwd = 2)
legend(x=230, y=3.8, legend = paste(round(strikes,4)), col = 1:3, lty=1, cex=0.8)
dev.off()

# For the whole smiles
rb <- roughbergomiclass(N = 200000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = 123)
simulate.roughbergomiclass  <- function(...) simulate_rb_antimixed(...)
price.roughbergomiclass     <- function(...) price_rb_mixed(...)

strikes <- loaddate(date = "2006-10-30")$volgrid$k
cat(round(strikes,2), sep = ", ")
x <- c(1:10 * 40, 6:16 * 80)
y <- sapply(x, function (n) {print(n); fastprice(rb, t = 0.25, k = strikes, n = n, pricetype = "impvol", seed = 123)} )

# save(strikes, x, y, file = "data/cht7_converge_n_whole_smile.RData")
# load(file = "data/cht7_converge_n_whole_smile.RData")

# Normal plot
pdf("Images/cht7_plot6.pdf")
par(mfrow=c(1,1), mar = c(4.5,4.5,1,0.5))
plot(x*0.25, y[1,] - mean(y[1,]), type = "l", ylim = c(-0.05,0.05))
for(i in 2:nrow(y)) {
     lines(x*0.25, y[i,] - mean(y[i,]), col = i)
}
dev.off()



# Adjusted plot
pdf("Images/cht7_plot7.pdf")
par(mfrow=c(1,1), mar = c(4.5,4.5,1,0.5))
plot(x*0.25, abs(y[1,] - y[1,ncol(y)]), type = "l", ylim = c(-0.001,0.05), xlab = "time steps (n/4)", ylab = "absolute error")
for(i in 2:nrow(y)) {
     lines(x*0.25, abs(y[i,] - y[i,ncol(y)]), col = i)
}
legend(x=250, y=0.05, legend = paste(round(strikes,2)), col = 1:nrow(y), lty=1, cex=0.8)
dev.off()

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

# save(rbclasses11, file = "data/cht7_volgrids_base.RData")

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

# save(rbclasses12, file = "data/cht7_volgrids_antithetic.RData")

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

# save(rbclasses13, file = "data/cht7_volgrids_conditional.RData")

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

# save(rbclasses14, file = "data/cht7_volgrids_controlled.RData")

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

# save(rbclasses15, file = "data/cht7_volgrids_mixed.RData")

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

# save(rbclasses16, file = "data/cht7_volgrids_antimixed.RData")

# Load
# load(file = "data/cht7_volgrids_base.RData")
# load(file = "data/cht7_volgrids_antithetic.RData")
# load(file = "data/cht7_volgrids_conditional.RData")
# load(file = "data/cht7_volgrids_controlled.RData")
# load(file = "data/cht7_volgrids_mixed.RData")
# load(file = "data/cht7_volgrids_antimixed.RData")

# Table
rbclasseslist <- list(rbclasses11,rbclasses12,rbclasses13,rbclasses14,rbclasses15,rbclasses16)
table <- gettablelist(rbclasseslist, exact = c(0.2961,0.2061,0.1576), names = c("base", "antithetic", "conditional", "controlled", "mixed", "antimixed"))
xtable(table, digits = 3)

# Plot
pdf("Images/cht7_plot3.pdf")
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

# save(exact2, file = "data/cht7_volgrids3_exact.RData")

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

# save(rbclasses21, file = "data/cht7_volgrids3_base.RData")

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

# save(rbclasses22, file = "data/cht7_volgrids3_antithetic.RData")

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

# save(rbclasses23, file = "data/cht7_volgrids3_conditional.RData")

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

# save(rbclasses24, file = "data/cht7_volgrids3_controlled.RData")

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

# save(rbclasses25, file = "data/cht7_volgrids3_mixed.RData")

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

# save(rbclasses26, file = "data/cht7_volgrids3_antimixed.RData")

# Load (t=0.25)
# load(file = "data/cht7_volgrids2_base.RData")
# load(file = "data/cht7_volgrids2_antithetic.RData")
# load(file = "data/cht7_volgrids2_conditional.RData")
# load(file = "data/cht7_volgrids2_controlled.RData")
# load(file = "data/cht7_volgrids2_mixed.RData")
# load(file = "data/cht7_volgrids2_antimixed.RData")
# load(file = "data/cht7_volgrids2_exact.RData")

# Load (t=1)
# load(file = "data/cht7_volgrids3_base.RData")
# load(file = "data/cht7_volgrids3_antithetic.RData")
# load(file = "data/cht7_volgrids3_conditional.RData")
# load(file = "data/cht7_volgrids3_controlled.RData")
# load(file = "data/cht7_volgrids3_mixed.RData")
# load(file = "data/cht7_volgrids3_antimixed.RData")
# load(file = "data/cht7_volgrids3_exact.RData")

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
xtable(t(table), digits = 3)

# plots
pdf("Images/cht7_plot4.pdf")
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses21, exact = exact2, xlim = c(-0.06,0.06))
dev.off()

pdf("Images/cht7_plot5.pdf")
par(mfrow = c(5,3), mar = c(2,2,1,0.5))
plotlist(rbclasses26, exact = exact2, xlim = c(-0.06,0.06))
dev.off()
