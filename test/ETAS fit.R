# Analyze Greenland (fit and bootstrap)
# Simulation study: a/b ratio approx. 1 and << 1
# Ogata '88 analysis
library(ETAS)
library(misc)
clean_up()


utsu = read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/UCPH/CIBS/Hawkes project/data/utsu.txt")
# tmp = as.POSIXct(paste0(utsu$YEAR,"/",utsu$MO,"/",utsu$DY," ", utsu$HR,":",utsu$MN), tz="UTC")
date = as.Date(paste0(utsu$YEAR,"/",utsu$MO,"/",utsu$DY))
time = chron::as.times(paste0(utsu$HR,":",utsu$MN,":00"))
mag  = utsu$MAG
shock = utsu$C # 0: main shock, 1: foreshoch, 2: aftershock

data("japan.quakes")
head(japan.quakes)
nrow(japan.quakes)

utsu = data.frame(date,time,long=NA, lat=NA, mag)
utsu$long = utsu$lat = runif(nrow(utsu),min = 0, max=1e-10)
head(utsu)
nrow(utsu)

catalog(utsu)
etas()