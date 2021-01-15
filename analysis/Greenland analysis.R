library(misc)
library(bootha)
misc::clean_up()
set.seed(123)

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
pp
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

# Plot fitted model
hi = 260
lo = 0
dt = 1/6 # hour resolution, i.e. 10 mins
x = seq(0,Tmax,dt)
y = intensity(x, pp, pars)
z = pp[pp > lo & pp < hi]

# h = hist(tw_hour, breaks=c(seq(0,260,dt)), ann=FALSE, border=NA, col=add.alpha('dodgerblue',.75))
h = hist(tw_hour, breaks=c(seq(0,260,dt)), plot=FALSE)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/UCPH/CIBS/Hawkes project/Greenland/fig/"

col1 = add.alpha('dodgerblue',.75)
col2 = add.alpha('tomato2',.85)
set_normal_plot_bg()



# pdf(paste0(plotlib,"greenland.pdf"), height=4, width=9)
    par(mfrow=c(1,1), mar=c(5,3,4,1), oma=c(0,0,0,0))
    plot(h$mids, h$counts, type='l', xlim=c(lo,hi), bty='n', col=col1, ann=FALSE, lwd=2)
    lines(x, y*dt, col=col2)
    mtext("# tweets/hour", side=2, line=2)
    mtext("hours since tweet",side = 1, line = 2.5)
    mtext("Greenland Twitter storm", side=3, line=2, cex=1.5)
    mtext("Between June 14-25th, 2019", side=3, line=.8, cex=1.15)
    legend("topright", c("Observed","Fitted Hawkes"), fill=c(col1,col2), bty='n')
# dev.off()

# points(z, rep(0,length(z)), pch=124)

L.exp = integrated_intensity(x, pp, pars)
s.exp = interpol(cbind(x,L.exp), pp)

s = s.exp
w = diff(s)
v = w/mean(w);mean(w);sd(w)


col1 = add.alpha('tomato2',.85)
col2 = add.alpha('steelblue',.5)

# pdf(paste0(plotlib,"diagnostics.pdf"), height=8, width=8)

    par(mfrow=c(2,2),mar=c(4,7,1,1), oma=c(0,0,3,0), bty='n')
    qqplot(w, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col1, pch=16, ann=FALSE); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2, las=0)
    mtext("Uncorrected rescaled waiting times w", 2, line=5.5, las=0)
    mtext("QQ plot", 3, line=2, las=0)
    ks_plot(w, col1 = col1, col2 = col2)
    mtext("KS plot", 3, line=2, las=0)
    
    # hist(w, breaks=seq(0,100,.01), prob=TRUE, col=col2, border=NA, ylim=c(0,1), ann=FALSE, xlim=c(0,1))
    # curve(dexp, add=TRUE, from=1e-8, col=col1, lwd=2)
    # mtext("Histogram of w", 3, line=2, las=0)
    
    
    qqplot(v, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col1, pch=16, ann=FALSE); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2, las=0)
    mtext("Mean corrected rescaled waiting times w", 2, line=5.5, las=0)
    ks_plot(v, col1 = col1, col2 = col2)
    
    # hist(v, breaks=seq(0,100,.01), prob=TRUE, col=col2, border=NA, ylim=c(0,1), ann=FALSE, xlim=c(0,1))
    # curve(dexp, add=TRUE, from=1e-8, col=col1, lwd=2)
    # mtext("Histogram of v", 3, line=2, las=0)

# dev.off()
ks.test(w,'pexp')
ks.test(v,'pexp')
pars

layout(1:3)
acf(diff(pp))
acf(w^2)
acf(v^2)


B    = 399
boot = bootstrap_hawkes(pp, pars, Tmax, B = B, type = 'fib', parametric = FALSE, ncores = 6)

boot.pars = matrix(nr=B, nc=length(pars))
for(i in 1:B){
    boot.pars[i,] = mle(boot[[i]], Tmax, mle.exp, dt=0)
    displayprogress(100*i/B)
}
boot.pars = as.data.frame(boot.pars)
names(boot.pars) = c("mu","n","B")
# boot.pars

f.mu = ecdf(boot.pars$mu)
f.n = ecdf(boot.pars$n)
f.B = ecdf(boot.pars$B)

f.mu(mle.exp[1])
f.n(mle.exp[2])
f.B(mle.exp[3])

layout(1:3)
plot(f.mu)
plot(f.n)
plot(f.B)

quantile(boot.pars$mu, c(.025,.975))
quantile(boot.pars$n, c(.025,.975))
quantile(boot.pars$B, c(.025,.975))

pars
x = seq(0,260,1)
layout(1)
plot(x, integrated_intensity(x, pp, pars), type='l')
for(i in 1:B){
  lines(x, integrated_intensity(x, boot[[i]], as.numeric(boot.pars[i,])), col=add.alpha('black',.25))
}
lines(x, integrated_intensity(x, pp, pars), type='l', col=add.alpha('red',.85))


