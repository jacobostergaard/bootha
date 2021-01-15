library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"


# Load the data
filelib   = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/criticality/"
fn        = paste0(filelib,"data/Greenland.Rda")
load(fn)

# To numeric (in seconds) and reset time to start at 0 (time of original tweet)
tw_time = as.numeric(Greenland$created_at)
tw_time = tw_time-min(tw_time)

# Find any timestamps within same second
tmp = as.data.frame(table(tw_time))
tmp = tmp[order(tmp$Freq,decreasing=TRUE),]


min(Greenland$created_at)

# Add a random (up to + .5s) time to timestamps within the same second to make obs unique
tw_time[tw_time %in% tmp$tw_time[tmp$Freq>1]] = tw_time[tw_time %in% tmp$tw_time[tmp$Freq>1]]+runif(sum(tmp$Freq[tmp$Freq>1]), min = 0,max = .5)
tw_time = sort(tw_time)

# Look at data in hourly timebins
range(tw_time/3600) # Timespan in hours of data
tw_hour = tw_time/3600

set_cobalt_plot_bg()
h = hist(tw_hour, breaks=c(seq(0,260,1)), ann=FALSE, border=NA, col=add.alpha('dodgerblue',.75))
mtext("# tweets/hour", side=2, line=2)
mtext("hours since tweet",side = 1, line = 2.5)
mtext("Greenland Twitter storm", side=3, line=2, cex=1.5)
mtext("Between June 14-25th, 2019", side=3, line=.8, cex=1.15)


range(diff(tw_time), na.rm = TRUE)

# Fit Hawkes

Tmax = 10*ceiling(max(tw_hour)/10)
pp = tw_hour[tw_hour<Tmax]

mu0 = length(pp)/Tmax
n0  = 1-(1/mu0)
B0  = 1

dt = 1
x = seq(0,Tmax,dt)

pars.exp = c(mu0,n0,B0)
pars.pow = c(mu0,n0,B0,1)

system.time({
  mle.exp = mle(pp, Tmax, pars.exp, dt=0) # fit exp kernel, this is instantaneous
})

datalib   = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
fn        = paste0(datalib,"mle_power_kernel.Rda")

system.time({
  if(file.exists(fn)){
    load(fn)
  } else{
    mle.pow = mle(pp, Tmax, pars.pow, dt=1) # fit pow kernel, this takes a little more than 1 hour...  
  }
})

mle.exp
mle.pow[2]/mle.pow[4]
mle.exp[1]
mle.pow[1]

pars = mle.exp



L.exp = integrated_intensity(x, pp, pars)
s.exp = interpol(cbind(x,L.exp), pp)

s = s.exp
w = diff(s)
v = w/mean(w);mean(w);sd(w)



set_normal_plot_bg()

layout(matrix(c(1,1,2,3,4,5), nr=3, nc=2, byrow=TRUE))
par(mar=c(3,3,1,1), oma=c(0,0,3,0))
plot_acf(diff(pp), lag=100)
mtext("ACF for Greenland")
# par(mfrow=c(2,2))

plot_acf(w, lag=100); mtext("Uncorrected waiting times", 3, line=0)
plot_acf(v, lag=100); mtext("Mean corrected waiting times", 3, line=0)
plot_acf(w^2, lag=100); mtext("Uncorrected waiting times squared", 3, line=0)
plot_acf(v^2, lag=100); mtext("Mean corrected waiting times squared", 3, line=0)
