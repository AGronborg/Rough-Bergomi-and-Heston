##############################$
##### IMPORTS & LIBRARIES #####
##############################$

### IMPORTS ###
source("thesislib.R")

### OPTIONS ###
par(mar=c(2,2,1,0.5)) # bot, left, top, right
par(mfrow=c(3,2))

####################$
##### PROCESSES #####
####################$

TT <- 3; n  <- 100; H <- 0.2

plot(sim_bm             (n = n, TT = TT),                                  type = "l")
plot(sim_fbm            (n = n, TT = TT, H = H),                           type = "l")
plot(sim_volterra       (n = n, TT = TT, H = H),                           type = "l")
plot(sim_volterra_hybrid(n = n, TT = TT, H = H, kappa = 0, b = "fwd"),     type = "l")
plot(sim_volterra_hybrid(n = n, TT = TT, H = H, kappa = 3, b = "optimal"), type = "l")

paths <- sim_bm_volterra (n = n, TT = TT, H = H, rho = -0.9)
plot(rbind(paths$bm,paths$volterra))

##################$
##### PRICING #####
##################$

### HESTON (LIPTON) FAST ###
priceheston(t = 1, k = -0.2, lambda = 0.2, vbar = 0.3^2, v0 = 0.2^2, eta = 0.02, rho = 0)

### HESTON (LIPTON) ###
hclass <- hestonclass(lambda = 0.2, vbar = 0.09, v0 = 0.04, eta = 0.02, rho = 0)
hclass <- setsimgrid(hclass, loaddate()$volgrid$t, loaddate()$volgrid$k)
hclass <- price_heston_closedform_lipton(hclass)

pricetime(hclass)
plot(hclass)

### ROUGH BERGOMI (MIXED) FAST ###
pricerb(t = 0.25, k = c(-0.1787,0.0000,0.1041), a = -0.43, rho = -0.9, eta = 1.9, xi = 0.235^2)

### ROUGH BERGOMI (MIXED) ###
rbclass <- roughbergomiclass(n = 1000, N = 1000, a = -0.43, rho = -0.90, eta = 1.9, xi = 0.235^2)
rbclass <- setsimgrid(rbclass, loaddate()$volgrid$t, loaddate()$volgrid$k)
rbclass <- simulate_rb_antimixed(rbclass)
rbclass <- price_rb_mixed(rbclass, estimator = "controlled")

simpricetime(rbclass)
plot(rbclass)
plot(rbclass$paths, variables = c("V","S1") , paths = c(1,501))

######################$
##### CALIBRATION #####
######################$

### DATA ###
dates <- loaddates_range(from = "2006-10-30", to = "2009-07-31", step = 300)
dates <- subset(dates, t = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3))

### HESTON (LIPTON) ###
hclass               <- hestonclass()
simulate.hestonclass <- function(...) identity(...)
price.hestonclass    <- function(...) price_heston_closedform_lipton(...)

hcotclass  <- calovertimeclass(dates, hclass)
hcotclass  <- calibrate_first(hcotclass, lambda = 1.15, vbar = 0.04, v0 = 0.04, eta = 0.39, rho = 0, control = list(trace = 1))
hcotclass  <- setvars_tofirst(hcotclass)
hcotclass  <- calibrate(hcotclass, v0 = hcotclass$simclasses[[1]]$vars$v0, eta = hcotclass$simclasses[[1]]$vars$eta, skipfirst = TRUE, control = list(trace = 1))
summary(hcotclass)
plot(hcotclass, ylim = c(0,0.6))

### ROUGH BERGOMI (MIXED) ###
rbclass                    <- roughbergomiclass(n = 500, N = 5000, seed = 123)
simulate.roughbergomiclass <- function(...) simulate_rb_antimixed(...)
price.roughbergomiclass    <- function(...) price_rb_mixed(...)
calibrate.roughbergomiclass<- function(...) calibrate_rb_mixed(...)

rbcotclass  <- calovertimeclass(dates, rbclass)
rbcotclass  <- calibrate_first(rbcotclass, a = -0.43, rho = -0.9, eta = 1.9, xi = 0.235^2, trackpars = TRUE, plottrace = TRUE, control = list(trace = 1, maxit = 20, REPORT = 1))
rbcotclass  <- setvars_tofirst(rbcotclass)
rbcotclass  <- calibrate(rbcotclass, eta = rbcotclass$simclasses[[1]]$vars$eta, xi = rbcotclass$simclasses[[1]]$vars$xi, skipfirst = TRUE, trackpars = TRUE, control = list(trace = 1, maxit = 20, REPORT = 1))
summary(rbcotclass)
plot(rbcotclass, ylim = c(0,0.6))

### PLOTS ###
multiplot_calovertime(list(hcotclass,rbcotclass), ylim = c(0,0.8))

### TABLES ###
xtable_h_rb(hcotclass, rbcotclass, c(1,1,2,2,1,1,1,2,2), digits = 3)