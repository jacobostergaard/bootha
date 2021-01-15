library(misc)
misc::clean_up()
set.seed(1234)

dji <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/DJI.csv")
dji$Date = as.Date(dji$Date)
dji$log.return = c(NA,diff(log(dji$Close)))

date.from = as.Date("1994-01-01")
date.to = as.Date("2010-12-31")

idx.dji = dji$Date>date.from & dji$Date < date.to
idx = c(1,5,7,8)
dji.tmp = dji[idx.dji,idx]
pfc = dji.tmp
names(pfc) = c("date", "close.dji", "volume.dji", "log.return.dji")

qts = quantile(pfc$log.return, probs = c(.01,.99), na.rm = TRUE)
qts = quantile(pfc$log.return, probs = c(.1,.9), na.rm = TRUE)

idx = which(pfc$log.return < qts[1] | pfc$log.return > qts[2])

x = as.numeric(pfc$date[idx])
x = x-min(x)+1
y = pfc$log.return[idx]

# Parameters from Embrechts
mu = 0.018
n = .74
d = 0.021
pp = as.numeric(x[y<0])
pars = c(mu,n,d)
Tmax = max(x)

dowjones = list(pp.neg = as.numeric(x[y<0]), pp.pos = as.numeric(x[y>0]), Tmax = Tmax, pars = pars)
dowjones 

datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
fn = paste0(datalib,"dowjones.Rda")
save(dowjones, file=fn)

