library(ppboot)
library(misc)
misc::clean_up()

dji <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/DJI.csv")
dji$Date  = as.Date(dji$Date)
dji$log.return = c(NA,diff(log(dji$Close)))
date.from = as.Date("1994-01-01")
date.to   = as.Date("2010-12-31")

idx = dji$Date>date.from & dji$Date < date.to
dji = dji[idx,c(1,5,7,8)]
names(dji) = c("date", "close.dji", "volume.dji", "log.return.dji")
head(dji)

qts = quantile(dji$log.return, probs = c(.1,.9), na.rm = TRUE)
idx = which(dji$log.return < qts[1] | dji$log.return > qts[2])

x = dji$date[idx]
y = dji$log.return[idx]
z = x[y<0]
yrs = unique(substr(z,1,4))
n.yrs = length(yrs)
# yrs = c(tmp,as.character(as.numeric(tmp[length(tmp)])+1)) # add a year to capture last period
# yrs = as.Date(paste0(yrs,"-01-01"))
pps = list()
for(i in 1:n.yrs){
  idx = substr(z,1,4) == yrs[i]
  pps[[i]] = as.numeric(z[idx]-as.Date(paste0(yrs[i],"-01-01")))
}
head(pps)

yrs = as.numeric(yrs)
par(mfrow=c(2,1), mar=c(3,2,1,1), oma=c(0,3,0,0), las=1, bty='n')
layout(1:2, heights = c(2,3))
plot(cbind(yrs,unlist(lapply(pps,length))), type='l', xaxt='n', yaxt='n')
abline(v=yrs, lty=3, col=add.alpha('black',.5))
abline(h=10*(1:6), lty=3, col=add.alpha('black',.5))
mtext("Events per year", side=2, line=3,las=0)
axis(side=1, at=yrs, labels=yrs, tick=FALSE, cex.axis=.75)
axis(side=2, at=10*(1:6), tick=FALSE, cex.axis=.75)


plot(0,0, type='l', xlim=c(0,365), ylim=c(n.yrs,1), xaxt='n', yaxt='n')
for(i in 1:length(yrs)){
  tmp = pps[[i]]
  points(tmp, rep(i, length(tmp)), pch=124)
}
mths = cumsum(c(0,31,28,31, 30, 31,30,31,30,31,31,30,31))
abline(v=mths, lty=3, col=add.alpha('black',.5))
lbls = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
axis(side=1, at=mths[-length(mths)]+diff(mths)/2, labels = lbls, tick=FALSE, cex.axis=.75)
axis(side=2, at=(n.yrs:1), rev(yrs), las=1, cex.axis=.75, tick=FALSE)
abline(h=1:(n.yrs+1)-.5, lty=3, col=add.alpha('black',.5))
mtext("Events each year", side=2, line=3, las=0)


