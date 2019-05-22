##############################$
##### IMPORTS & LIBRARIES #####
##############################$

### IMPORTS ###
source("thesislib.R")
library(beepr)

####################$
##### CHAPTER 1 #####
####################$

##########################$
##### CHAPTER 2- DATA #####
##########################$

# Data
data  <- loaddate(date = "2006-10-30")

length(data$volgrid$k)
length(data$volgrid$t)
13*14
# Strikes
table           <- rbind(paste(exp(data$volgrid$k)*100, "%", sep = ""), round(data$volgrid$k,2))
rownames(table) <- c("K", "k")

xtable(table, align = "c|ccccccccccccc")

# Maturities
table           <- rbind(paste(round(data$volgrid$t*12), "M", sep = ""), data$volgrid$t)
rownames(table) <- c("Months", "t")

xtable(table, align = paste("c|cccccccccccccc")) 

# Volatility surface
table           <- data.frame(data$volgrid$impvol)
rownames(table) <- paste(exp(data$volgrid$k)*100, "%", sep = "")
colnames(table) <- data$tenor
                       
xtable(table, align = "c|cccccccccccccc", digits = 3)

#############################$
##### CHAPTER 3 - HESTON #####
#############################$

# Pricing times
hclass         <- hestonclass(n = 100, N = 10000, TT = 1, lambda = 1.3253, vbar = 0.0354, v0 = 0.0174, eta = 0.3877, rho = -0.7165)
ptimegatheral  <- benchmodel(hclass, t = 1, k = 0, times = 100, unit = "millisecs", simfunc = identity, pricefunc = price_heston_closedform_gatheral)
ptimelipton    <- benchmodel(hclass, t = 1, k = 0, times = 100, unit = "millisecs", simfunc = identity, pricefunc = price_heston_closedform_lipton)
ptimemc        <- benchmodel(hclass, t = 1, k = 0, times = 100, unit = "millisecs", simfunc = simulate_heston, pricefunc = price_standard)


# Prices
hclass     <- hestonclass(n = 100, N = 10000, TT = 1, lambda = 1.3253, vbar = 0.0354, v0 = 0.0174, eta = 0.3877, rho = -0.7165)
pgatheral  <- replicate(fastprice(hclass, t = 1, k = 0, simfunc = identity, pricefunc = price_heston_closedform_gatheral), n = 100)
plipton    <- replicate(fastprice(hclass, t = 1, k = 0, simfunc = identity, pricefunc = price_heston_closedform_lipton), n = 100)
pmc        <- replicate(fastprice(hclass, t = 1, k = 0, simfunc = simulate_heston, pricefunc = price_standard), n = 100)

# Table 1
table <- data.frame(Runtime = c(ptimegatheral, ptimelipton, ptimemc),
                    Estimate = c(mean(pgatheral),mean(plipton),mean(pmc)),
                    Error = c(rse(pgatheral), rse(plipton), rse(pmc))) 
colnames(table) <- c("Run time", "Estimate", "RSE")
rownames(table) <- c("Gatheral", "Lipton", "Monte Carlo")
xtable(table, digits = 3)

#################################################$
##### CHAPTER 4 - FRACTIONAL BROWNIAN MOTION #####
#################################################$

# fbm plot 1 (different H)
H <- c(0.1,0.3,0.5,0.7,0.9,0.95)
paths <- lapply(H, sim_fbm, n = 500, TT = 1, N = 1, simseed = 123)

pdf("../Thesis LaTeX/Images/cht4_plot1.pdf")
par(mfrow = c(3,2), mar = c(2.5,2.5,1,1))
lapply(paths, plot, x = 0:500/500, type = "l", xlab = "", ylab = "")
dev.off()

# fbm plot 2 (scaling)
path1 <- sim_fbm(H = 0.3, n = 500, TT = 1,  N = 1, simseed = 123)
path2 <- sim_fbm(H = 0.3, n = 50, TT = 10, N = 1, simseed = 123)

pdf("../Thesis LaTeX/Images/cht4_plot2.pdf")
par(mfrow = c(1,1), mar = c(2.5,2.5,1,1))
plot(seq(0,10,length.out = 501), path1*10^0.3, type = "l", lwd = 3)
lines(seq(0,10,length.out = 501), path2, col = 2, lwd = 1)
dev.off()

# Correlation of increaments
paths <- sim_fbm(H = 0.3, n = 250, TT = 2, N = 50000, simseed = 123)
incre <- apply(paths, 1, function(x) c(x[251]-x[1],x[501]-x[251])) 
cor(incre[1,],incre[2,])

4^(0.3-0.5)-1

# Distribution of increments
paths <- sim_fbm(H = 0.3, n = 250, TT = 2, N = 50000, simseed = 123)
incre1 <- apply(paths, 1, function(x) x[251]-x[1])
incre2 <- apply(paths, 1, function(x) x[501]-x[251])

x <- seq(-3,3, length.out = 100)
pdf("../Thesis LaTeX/Images/cht4_plot3.pdf")
plot(x, dnorm(x), type = "l", ylab = "density")
lines(density(incre1)$x, density(incre1)$y, col = 2)
lines(density(incre2)$x, density(incre2)$y, col = 3)
dev.off()

####################################$
##### CHAPTER 5 - ROUGH BERGOMI #####
####################################$

### DATA ###
empirical <- loaddate(date = "2006-10-30")
k <- seq(-0.5,0.5, by = 0.05)

# ##### EXACT SIMULATION #####
# # t = 0.041
# rbclass1 <- roughbergomiclass(n = round(50/0.041), N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2)
# rbclass1 <- setsimgrid(rbclass1, t = 0.041, k = c(log(0.6),log(0.65),empirical$volgrid$k))
# rbclass1 <- simulate_rb_exact(rbclass1)
# rbclass1 <- price_standard(rbclass1)
# grid1    <- rbclass1$simgrid
# 
# # t = 1
# rbclass2 <- roughbergomiclass(n = 50, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2)
# rbclass2 <- setsimgrid(rbclass2, t = 1, k = c(log(0.6),log(0.65),empirical$volgrid$k, log(1.4), log(1.5), log(1.6)))
# rbclass2 <- simulate_rb_exact(rbclass2)
# rbclass2 <- price_standard(rbclass2)
# grid2    <- rbclass2$simgrid

##### EXACT SIMULATION (SPLITTING)#####
# t = 0.25
rbclass1 <- roughbergomiclass(n = 2000, N = 20000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
rbclass1 <- setsimgrid(rbclass1, t = 0.25, k = k)
grid1    <- splitprice(rbclass1, times = 50, sim_func = function(x) simulate_rb_exact(x))

save(grid1, file = "data/cht5_exact1.RData")

# t = 1
rbclass2 <- roughbergomiclass(n = 500, N = 20000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
rbclass2 <- setsimgrid(rbclass2, t = 1, k = k)
grid2    <- splitprice(rbclass2, times = 50, sim_func = function(x) simulate_rb_exact(x))

save(grid2, file = "data/cht5_exact2.RData")

##### PLOTS #####
pdf("../Thesis LaTeX/Images/cht5_plot1.pdf")
plot(grid1, type = "l", mar = c(4.5,4.5,2,2))
dev.off()

pdf("../Thesis LaTeX/Images/cht5_plot2.pdf")
plot(grid2, type = "l", mar = c(4.5,4.5,2,2))
dev.off()

####################################$
##### CHAPTER 6 - HYBRID SCHEME #####
####################################$

### DATA ###
empirical <- loaddate(date = "2006-10-30")
k <- seq(-0.5,0.5, by = 0.05)

##### PRICING #####
# # t = 0.041
# rbclass   <- roughbergomiclass(n = round(500/0.041), N = 1000000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
# rbclass   <- setsimgrid(rbclass, t = 0.041, k = c(log(0.6),log(0.65),empirical$volgrid$k))
# 
# price.roughbergomiclass <- function(x) price_standard(x)
# grid11 <- fastprice(rbclass, simfunc = function(x) simulate_rb_exact(x))
# grid12 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 0, b = "fwd"))
# grid13 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 0, b = "optimal"))
# grid14 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 1, b = "optimal"))
# grid15 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 3, b = "optimal"))
# 
# # t = 1
# rbclass   <- roughbergomiclass(n = 500, N = 1000000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
# rbclass   <- setsimgrid(rbclass, t = 1, k = c(log(0.6),log(0.65),empirical$volgrid$k, log(1.4), log(1.5), log(1.6)))
# 
# price.roughbergomiclass <- function(x) price_standard(x)
# grid21 <- fastprice(rbclass, simfunc = function(x) simulate_rb_exact(x))
# grid22 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 0, b = "fwd"))
# grid23 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 0, b = "optimal"))
# grid24 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 1, b = "optimal"))
# grid25 <- fastprice(rbclass, simfunc = function(x) simulate_rb_standard(x, kappa = 3, b = "optimal"))

##### PRICING (SPLITTING) #####
# t = 0.25
rbclass   <- roughbergomiclass(n = 500, N = 20000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
rbclass   <- setsimgrid(rbclass, t = 0.25, k = k)

price.roughbergomiclass <- function(x) price_standard(x)
grid11     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_exact(x))
grid12     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 0, b = "fwd"))
grid13     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 0, b = "optimal"))
grid14     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 1, b = "optimal"))
grid15     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 3, b = "optimal"))

save(grid11, file = "data/cht6_hybrid11.RData")
save(grid12, file = "data/cht6_hybrid12.RData")
save(grid13, file = "data/cht6_hybrid13.RData")
save(grid14, file = "data/cht6_hybrid14.RData")
save(grid15, file = "data/cht6_hybrid15.RData")

# t = 1
rbclass   <- roughbergomiclass(n = 500, N = 20000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = -1)
rbclass   <- setsimgrid(rbclass, t = 1, k = k)

price.roughbergomiclass <- function(x) price_standard(x)
grid21     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_exact(x))
grid22     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 0, b = "fwd"))
grid23     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 0, b = "optimal"))
grid24     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 1, b = "optimal"))
grid25     <- splitprice(rbclass, times = 50, sim_func = function(x) simulate_rb_standard(x, kappa = 3, b = "optimal"))

save(grid21, file = "data/cht6_hybrid21.RData")
save(grid22, file = "data/cht6_hybrid22.RData")
save(grid23, file = "data/cht6_hybrid23.RData")
save(grid24, file = "data/cht6_hybrid24.RData")
save(grid25, file = "data/cht6_hybrid25.RData")

##### TABLES #####
gridlist  <- list(grid11, grid12, grid13, grid14, grid15)
errors    <- sapply(gridlist, function(x) abs(x$impvol - grid11$impvol))
relerrors <- sapply(gridlist, function(x) abs(x$impvol - grid11$impvol)/grid11$impvol)

colnames(relerrors) <- c("exact", "fwd, kappa = 0", "optimal, kappa = 0", "optimal, kappa = 1", "optimal, kappa = 3")
rownames(relerrors) <- round(k, 2)
xtable(relerrors*100, digits = 3)

gridlist  <- list(grid21, grid22, grid23, grid24, grid25)
errors    <- sapply(gridlist, function(x) abs(x$impvol - grid21$impvol))
relerrors <- sapply(gridlist, function(x) abs(x$impvol - grid21$impvol)/grid21$impvol)

colnames(relerrors) <- c("exact", "fwd, kappa = 0", "optimal, kappa = 0", "optimal, kappa = 1", "optimal, kappa = 3")
rownames(relerrors) <- k
xtable(relerrors*100, digits = 3)

##### PLOTS #####
pdf("../Thesis LaTeX/Images/cht6_plot1.pdf")
par(mfrow = c(2,1), mar = c(4.5,4.5,2,2))
plotmultiplevolgrids(list(grid11, grid12, grid13, grid14, grid15), pricetype = "impvol", mfrow = FALSE, y_min = 0, lwd = 2)
plotmultiplevolgrids(list(grid21, grid22, grid23, grid24, grid25), pricetype = "impvol", mfrow = FALSE, y_min = 0, lwd = 2)
dev.off()

pdf("../Thesis LaTeX/Images/cht6_plot2.pdf")
par(mfrow = c(2,1), mar = c(4.5,4.5,2,2))
plotmultiplevolgrids(list(grid11, grid12, grid13, grid14, grid15), pricetype = "prices", mfrow = FALSE, y_min = 0)
plotmultiplevolgrids(list(grid21, grid22, grid23, grid24, grid25), pricetype = "prices", mfrow = FALSE, y_min = 0)
dev.off()

##### PRICING TIME #####

rbclass <- roughbergomiclass(a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, N = 10)
n       <- 1:10*10

time1 <- 1/10*sapply(n, function(n) benchmodel(rbclass, n = n, times = 1, simfunc = function(x) simulate_rb_exact(x)))
time2 <- 1/10*sapply(n, function(n) benchmodel(rbclass, n = n, times = 1, simfunc = function(x) simulate_rb_standard(x, kappa = 0, b = "fwd")))
time3 <- 1/10*sapply(n, function(n) benchmodel(rbclass, n = n, times = 1, simfunc = function(x) simulate_rb_standard(x, kappa = 0, b = "optimal")))
time4 <- 1/10*sapply(n, function(n) benchmodel(rbclass, n = n, times = 1, simfunc = function(x) simulate_rb_standard(x, kappa = 1, b = "optimal")))
time5 <- 1/10*sapply(n, function(n) benchmodel(rbclass, n = n, times = 1, simfunc = function(x) simulate_rb_standard(x, kappa = 3, b = "optimal")))

par(mfrow=c(1,1))
plot(rbind(time1,time2,time3,time4,time5), xvalues = 1:10*4, paths = 1:5, xlab = "n", ylab = "Time per path")

# Regressions
y <- time1
x <- n^3
reg1 <- lm(y ~ x)
plot(x = n, y = y)
lines(x = n, y = fitted(reg1))

y <- time5
x <- n*log(n)
reg2 <- lm(y ~ x)
plot(x = n, y = y)
lines(x = n, y = fitted(reg2))

pred1 <- predict(reg1, list(x = 312^3))
pred2 <- predict(reg1, list(x = 312*log(312)))

pred1/pred2

rbclass <- roughbergomiclass(n = 312, N = 10, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2)
rbclass <- setsimgrid(rbclass, t = 1, k = 0)
rbclass <- simulate_rb_exact(rbclass)
rbclass <- price_standard(rbclass)
sim1 <- simpricetime(rbclass)
rbclass <- simulate_rb_standard(rbclass, kappa = 3, b = "optimal")
rbclass <- price_standard(rbclass)
sim2 <- simpricetime(rbclass)

as.numeric(sim1)/as.numeric(sim2)

# Tables
table <- rbind(time1,time2,time3,time4,time5)
rownames(table) <- c("exact", "fwd, kappa = 0", "optimal, kappa = 0", "optimal, kappa = 1", "optimal, kappa = 3")
colnames(table) <- n

xtable(table*10^3, align = "r|rrrrrrrrrr", digits = 3)

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
load(file = "data/cht7_converge_n.RData")

# Absolute error
conf <- 1/sqrt(200) * c(0.0055,0.0027,0.0026) * qnorm(0.975) # 95% confidence interval (error)

pdf("../Thesis LaTeX/Images/cht7_plot1.pdf")
par(mfrow=c(1,1))
plot( x/4, abs(y[1,]-0.2961), type = "l", col = 1, ylim = c(0,0.02), xlab = "time steps", ylab = "absolute error")
lines(x/4, abs(y[2,]-0.2061), type = "l", col = 2)
lines(x/4, abs(y[3,]-0.1576), type = "l", col = 3)
abline(h = conf[1:3], col = 1:3, lwd = 2)
dev.off()

# Relative error
conf <- 1/sqrt(200) * c(0.0055,0.0027,0.0026) * qnorm(0.975) / c(0.2961,0.2061,0.1576) * 100 # 95% confidence interval (error)

pdf("../Thesis LaTeX/Images/cht7_plot2.pdf")
par(mfrow=c(1,1))
plot( x/4, abs(y[1,]-0.2961)/0.2961*100, type = "l", col = 1, ylim = c(0,4), xlab = "time steps", ylab = "relative absolute error")
lines(x/4, abs(y[2,]-0.2061)/0.2061*100, type = "l", col = 2)
lines(x/4, abs(y[3,]-0.1576)/0.1576*100, type = "l", col = 3)
abline(h = conf[1:3], col = 1:3, lwd = 2)
dev.off()

# For the whole smiles
rb <- roughbergomiclass(N = 200000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, seed = 123)
simulate.roughbergomiclass  <- function(...) simulate_rb_antimixed(...)
price.roughbergomiclass     <- function(...) price_rb_mixed(...)

strikes <- empirical$volgrid$k
x <- c(1:10 * 40, 6:16 * 80)
y <- sapply(x, function (n) {print(n); fastprice(rb, t = 0.25, k = strikes, n = n, pricetype = "impvol", seed = 123)} )

#save(strikes, x, y, file = "data/cht7_converge_n_whole_smile.RData")
load(file = "data/cht7_converge_n_whole_smile.RData")

# Normal plot
plot(x*0.25, y[1,] - mean(y[1,]), type = "l", ylim = c(-0.05,0.05))
for(i in 2:nrow(y)) {
     lines(x*0.25, y[i,] - mean(y[i,]), col = i)
}

# Adjusted plot
plot(x*0.25, y[1,] - mean(y[1,]), type = "l", ylim = c(-0.005,0.05))
for(i in 2:nrow(y)) {
     lines(x*0.25, y[i,] - y[i,ncol(y)], col = i)
}

##################################$
##### CHAPTER 8 - CALIBRATION #####
##################################$

##### HESTON #####
data                 <- loaddate(date = "2009-07-31")
data                 <- subset(data, t = data$volgrid$t[-c(1,2)])
heston               <- hestonclass()
simulate.hestonclass    <- function(...) identity(...)
price.hestonclass       <- function(...) price_heston_closedform_lipton(...)

heston <- setempgrid(heston, data)
heston <- calibrate(heston, lambda = 1.15, vbar = 0.04, v0 = 0.04, eta = 0.39, rho = 0)

# Table
getvars(heston, digits = 3, calinfo = TRUE)
xtable(getvars(heston, fortable = TRUE, calinfo = TRUE), digits = 3)

# Plot
pdf("../Thesis LaTeX/Images/cht8_plot1.pdf")
par(mar=c(2.3,2,1,1))
plot(heston, ylim = c(0,0.4)) 
dev.off()

##### ROUGH BERGOMI #####
data                       <- loaddate(date = "2009-07-31")
data                       <- subset(data, t = data$volgrid$t[-c(1,2)])
rb                         <- roughbergomiclass(n = 312*4, N = 1000, seed = 123)
simulate.roughbergomiclass  <- function(...) simulate_rb_antimixed(...)
price.roughbergomiclass     <- function(...) price_rb_mixed(...)

rb   <- setempgrid(rb, data)
rb   <- calibrate_rb_mixed(rb, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, control = list(maxit = 20, trace = 1, REPORT = 1), trackpars = TRUE)

# Table
getvars(rb, digits = 3, calinfo = TRUE)
xtable( getvars(rb, fortable = TRUE, calinfo = TRUE) , digits = 3)

# Plot
pdf("../Thesis LaTeX/Images/cht8_plot2.pdf")
par(mar=c(2.3,2,1,1))
plot(rb) 
dev.off()

##### COMPARISON #####

# Plot
pdf("../Thesis LaTeX/Images/cht8_plot3.pdf")
par(mar=c(2.3,2,1,1))
plotmultiplevolgrids(list(class0 = data$volgrid, class1 = heston$simgrid, class2 = rb$simgrid)) 
dev.off()

##########################################$
##### CHAPTER 9 - CALIBRATE OVER TIME #####
##########################################$

### DATA ###
dates <- loaddates_range(from = "2006-10-30", to = "2009-07-31", step = 300)
dates <- subset(dates, t = c(0.25,0.5,1,2))
showdates(dates)

### HESTON ###

# Calibration
hestoncls  <- hestonclass()
simulate.hestonclass    <- function(...) identity(...)
price.hestonclass       <- function(...) price_heston_closedform_lipton(...)
calibrate.hestonclass   <- function(...) calibrate.simulateclass(...)

hestoncot  <- calovertimeclass(dates, hestoncls)
hestoncot  <- calibrate_first(hestoncot, lambda = 1.15, vbar = 0.04, v0 = 0.04, eta = 0.39, rho = 0)
hestoncot  <- setvars_tofirst(hestoncot)
hestoncot  <- calibrate(hestoncot, v0 = hestoncot$simclasses[[1]]$vars$v0, eta = hestoncot$simclasses[[1]]$vars$eta, skipfirst = TRUE)

# Tables
getvars(hestoncot)
xtable(getvars(hestoncot))

summary(hestoncot)
xtable(summary(hestoncot))

# Plots
pdf("../Thesis LaTeX/Images/cht9_plot1.pdf")
par(mar=c(2.3,2,1,1))
plot(hestoncot)
dev.off()

### ROUGH BERGOMI ###

# Calibration
rbcls  <- roughbergomiclass(n = 312*4, N = 1000, seed = 123)
simulate.roughbergomiclass  <- function(...) simulate_rb_antimixed(...)
price.roughbergomiclass     <- function(...) price_rb_mixed(...)
calibrate.roughbergomiclass <- function(...) calibrate_rb_mixed(...)

rbcot  <- calovertimeclass(dates, rbcls)
rbcot  <- calibrate_first(rbcot, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2, control = list(maxit = 20, trace = 1, REPORT = 1), trackpars = TRUE)
rbcot  <- setvars_tofirst(rbcot)
rbcot  <- calibrate(rbcot, xi = rbcot$simclasses[[1]]$vars$xi, eta = rbcot$simclasses[[1]]$vars$eta, skipfirst = TRUE, control = list(maxit = 20, trace = 1, REPORT = 1), trackpars = TRUE)

# Tables
getvars(rbcot)
xtable(getvars(rbcot))

summary(rbcot)
xtable(summary(rbcot))

# Plots
pdf("../Thesis LaTeX/Images/cht9_plot2.pdf")
par(mar=c(2.3,2,1,1))
plot(rbcot)
dev.off()
