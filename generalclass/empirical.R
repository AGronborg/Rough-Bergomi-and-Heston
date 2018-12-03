###########################
### IMPORTS & LIBRARIES ###
###########################

##############
### IMPORT ###
##############

# R date format: yyyy-mm-dd 2009-09-30
# file date format: dd/mmm/yy 30/Jul/09
# Date from 2006-10-30 to 2009-07-31

finddaterow <- function(dataf, date) {
     match(format(date, "%d/%b/%y"),dataf[,2])
}

stringtonum <- function(string) {
     as.numeric(gsub(",", ".", string))
}

dateclass <- function(dataf, row = NULL, date = NULL) {
     d         <- list()
     attr(d, "class") <- "dateclass"
     
     if (!is.null(date)) row <- finddaterow(dataf,date) # Start row
     if (is.na(row)) stop("could not find date")
     
     nk        <- match("",dataf[(row+2),4:ncol(dataf)])-1  # number of strikes
     nt        <- match("",dataf[(row+3):nrow(dataf),2])-1  # number of maturities
     
     d$date      <- as.Date(dataf[row,2], "%d/%b/%y")
     d$tenor     <- dataf[(row+3):(row+3+nt-1),1]
     d$maturity  <- dataf[(row+3):(row+3+nt-1),2]
     d$S         <- stringtonum( dataf[row+1,2]                    )
     d$K         <- stringtonum( rev(dataf[row+1,4:(4+nk-1)])      )
     
     t         <- stringtonum( dataf[(row+3):(row+3+nt-1),3]       )
     k         <- log( rev(stringtonum(dataf[row+2,4:(4+nk-1)]))   ) # Alternatively: log( K/S )
     impvol    <- apply( dataf[(row+3):(row+3+nt-1),(4+nk-1):4], 2, stringtonum)
     
     d$volgrid        <- volatilitygrid(t = t, k = k, impvol = t(impvol))
     d$volgrid$prices <- getprices(d$volgrid
                                   )
     d$zerorate  <- stringtonum( dataf[(row+3):(row+3+nt-1),4+nk]  )
     d$divyield  <- stringtonum( dataf[(row+3):(row+3+nt-1),5+nk]  )
     d$forward   <- stringtonum( dataf[(row+3):(row+3+nt-1),6+nk]  )
     d$atmfvol   <- stringtonum( dataf[(row+3):(row+3+nt-1),7+nk]  )
     
     return(d)
}

subset.dateclass <- function(dclass, t = NULL, k = NULL) {     
     i               <- index(t, dclass$volgrid$t)
     
     dclass$volgrid  <- subset(dclass$volgrid, t = t, k = k)
     dclass$tenor    <- dclass$tenor[i]
     dclass$maturity <- dclass$maturity[i]
     dclass$K        <- dclass$K[i]
     dclass$zerorate <- dclass$zerorate[i]
     dclass$divyield <- dclass$divyield[i]
     dclass$forward  <- dclass$forward[i]
     dclass$atmfvol  <- dclass$atmfvol[i]
     
     return(dclass)
}

loaddate <- function(date = "2006-10-30") {
     dataf <- read.table("generalclass/volsurf_06_09.csv", sep = ";", as.is = TRUE)
     date  <- as.Date(date)

     row       <- finddaterow(dataf,date)                   # Start row
     if (is.na(row)) stop(paste("could not find date", date))
     
     dateclass(dataf, row)
}

loaddates <- function(dates) {
     dclass <- list()
     attr(dclass, "class") <- "datesclass"
     
     for (i in 1:length(dates)) dclass[[i]] <- loaddate(dates[i]) 
     dclass$info$nrobs <- length(dates)
     
     return(dclass)
}

loaddates_range <- function(from = "2006-10-30", to = "2009-07-30", step = 1, useoffset = TRUE) {
     dataf <- read.table("generalclass/volsurf_06_09.csv", sep = ";", as.is = TRUE)
     from <- as.Date(from)
     to   <- as.Date(to)
     
     datesdata <- list()
     attr(datesdata, "class") <- "datesclass"
     
     i <- 0
     while(from <= to) {
          row    <- finddaterow(dataf, from) # Start row
          offset <- 0
          while (is.na(row) && from <= to) {
               from   <- from + 1
               offset <- offset + 1
               row    <- finddaterow(dataf, from)
          }
          if (!is.na(row)) {
               i              <- i + 1
               datesdata[[i]] <- dateclass(dataf, row)
               if (useoffset && step > offset) from <- from + step - offset
               else if (useoffset) from <- from + 1
               else from <- from + step
          }
     }
     datesdata$info$nrobs <- i
     
     return(datesdata)
}

onlydates <- function(dates) {
     dates$info <- NULL
     return(dates)
}

showdates <- function(dates) {
     sapply(onlydates(dates), function(x) as.character(x$date))
}

subset.datesclass <- function(dsclass, t = NULL, k = NULL) {
     info           <- dsclass$info
     dsclass        <- onlydates(dsclass)
     dsclass        <- lapply(dsclass, function(x) subset(x, t, k))
     class(dsclass) <- "datesclass"
     dsclass$info   <- info
     return(dsclass)
}

plot.dateclass <- function(obj, style = c("normal","3d","extra"), ...) {
     style <- match.arg(style)
     
     if (style == "normal") {
          par(mfrow = c(4,4))
          plot(obj$volgrid, ...)
     } else if (style == "3d") {
          par(mfrow = c(1,1))
          plot(obj$volgrid, plot3d = TRUE, ...)
     } else if (style == "extra") {
          par(mfrow = c(2,2))
          plot(x = obj$volgrid$t, y = obj$zerorate, main = "zero rate", xlab = "t", ...)
          plot(x = obj$volgrid$t, y = obj$divyield, main = "dividend yield", xlab = "t", ...)
          plot(x = obj$volgrid$t, y = obj$forward, main = "forward", xlab = "t", ...)
          plot(x = obj$volgrid$t, y = obj$atmfvol, main = "ATM forward volatility", xlab = "t", ...)
     }
}