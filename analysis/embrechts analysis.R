# Embrechts
library(ppboot)
library(misc)
misc::clean_up()

gk <- function(x,a,b,l){
 
  tmp1 = l/(a*l+b)
  tmp2 = a+b*abs(x)
  out = tmp1*tmp2
  
  return(out)
}

wj <- function(x,d){
  out = d*exp(-d*x)
  return(out)
}

fk <- function(x,l){
  out = dexp(x,l)
  return(out)
}

cond.int1 <- function(x,N,M,n,V,d,a,b,l){
  # x=9000; m=1; i=1; j=1

  k   = length(n)   # number of Hawkes components
  out = numeric(k)  # output at time x is a vector of length k
  
  for(m in 1:k){
    out.tmp = 0
    for(i in 1:k){
      Nk  = N[[i]]
      Mk  = M[[i]]
      
      Nt  = sum(Nk < x)
      tmp = 0
      if(Nt>0){
        for(j in 1:Nt){
          w.tmp = wj(x-Nk[j],d[i])
          g.tmp = gk(Mk[j], a[i], b[i], l[i])
          tmp   = tmp+w.tmp*g.tmp
        }  
        v.tmp = V[m,i]*tmp
      } else{
        v.tmp = 0
      }
      out.tmp = out.tmp + v.tmp
    }  
    out[m] = n[m]+out.tmp
  }
  
  return(out)
}

cond.int <- function(x,N,M,n,V,d,a,b,l){
  out = matrix(nr=length(x), nc=length(n))
  for(i in 1:length(x)){
    out[i,] = cond.int1(x[i],N,M,n,V,d,a,b,l)
  }
  return(out)
}


dji <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/DJI.csv")
ndx <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/NDX.csv")
spc <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/SPC.csv")


dji$Date = as.Date(dji$Date)
ndx$Date = as.Date(ndx$Date)
spc$Date = as.Date(spc$Date)

range(dji$Date)
range(ndx$Date)
range(spc$Date)

dji.col = add.alpha('steelblue',.75)
ndx.col = add.alpha('tomato2',.75)
spc.col = add.alpha('darkgoldenrod2',.75)

spc$log.return = c(NA,diff(log(spc$Close)))
ndx$log.return = c(NA,diff(log(ndx$Close)))
dji$log.return = c(NA,diff(log(dji$Close)))

date.from = as.Date("1994-01-01")
date.to = as.Date("2010-12-31")

idx.dji = dji$Date>date.from & dji$Date < date.to
idx.ndx = ndx$Date>date.from & ndx$Date < date.to
idx.spc = spc$Date>date.from & spc$Date < date.to

idx = c(1,5,7,8)
dji.tmp = dji[idx.dji,idx]
ndx.tmp = ndx[idx.ndx,idx]
spc.tmp = spc[idx.spc,idx]

pfc = dji.tmp
names(pfc) = c("date", "close.dji", "volume.dji", "log.return.dji")

head(pfc)

plot(pfc$date, pfc$log.return, type='l')

qts = quantile(pfc$log.return, probs = c(.01,.99), na.rm = TRUE)
qts = quantile(pfc$log.return, probs = c(.1,.9), na.rm = TRUE)
abline(h=qts, lty=3)

idx = which(pfc$log.return < qts[1] | pfc$log.return > qts[2])

x = as.numeric(pfc$date[idx])
x = x-min(x)
y = pfc$log.return[idx]

par(bty='n', mar=c(2,2,2,2), oma=c(0,0,0,0), mfrow=c(2,1))

plot(x[y<0], y[y<0], type='h')
plot(x[y>0], y[y>0], type='h')


N = list(N1=as.numeric(x[y<0]), N2=as.numeric(x[y>0]))
M = list(M1=as.numeric(y[y<0]), M2=as.numeric(y[y>0]))

V = matrix(c(.74, 0.83, 0,0), nr=2, nc=2)
max(abs(eigen(V)$values)) # spectral radius
n = c(0.018,0.012)
a = c(1,1)
b = c(47,74)
l = c(109,122)
d = c(0.021,0.021)

dt = 1
x   = 1:max(as.numeric(x))
x   = seq(0,6144,dt)
system.time({
  lam = cond.int(x,N,M,n,V,d,a=c(1,1),b,l)  
})
Lam = apply(lam,2,cumsum)*dt
L = data.frame(x=x,y1=Lam[,1],y2 = Lam[,2])

s1 = interpol(as.matrix(L[,c(1,2)]), x=N[[1]])
s2 = interpol(as.matrix(L[,c(1,3)]), x=N[[2]])
w1 = diff(s1)
w2 = diff(s2)
v1 = w1/mean(w1);mean(w1)
v2 = w2/mean(w2);mean(w2)

sd(w1)
sd(w2)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/UCPH/CIBS/Hawkes project/Embrechts/fig/"

pdf(paste0(plotlib,"embrechts1.pdf"), height=5, width=7)
    par(mfrow=c(4,1), bty='n', oma=c(1,0,0,0), mar=c(2,5,1,1), las=1)
    layout(1:4,heights=c(8,8,2,2))
    tmp = as.Date(paste0(seq(1995,2010,5),"-01-01"))
    
    col1 = add.alpha('steelblue',.75)
    col2 = add.alpha('tomato2',.75)
    plot(N[[2]],M[[2]], type='h', col=col1, ylim=c(-.1,.1), xaxt='n', ann=FALSE)
    lines(N[[1]],M[[1]], type='h', col=col2)
    abline(h=seq(-.1,.1,.025), col=add.alpha('gray',.5))
    abline(h=0)
    axis(side = 1, at=as.numeric(tmp)-min(as.numeric(pfc$date[idx])), labels=substr(tmp,1,4))
    legend("topright",c("j=1","j=2"), col=c(col2,col1), bty='n', lty=1, lwd=2)
    mtext("Excess log-returns", side=2, line=3.5, las=0, cex=.75)
    
    plot(x,lam[,2], type='l', col=col1, xaxt='n', ylim=c(0,.5), ann=FALSE)
    lines(x,lam[,1], col=col2)
    axis(side = 1, at=as.numeric(tmp)-min(as.numeric(pfc$date[idx])), labels=substr(tmp,1,4))
    abline(h=seq(0,.5,.1), col=add.alpha('gray',.5))
    mtext("Estimated intensity", side=2, line=3.5, las=0, cex=.75)
    
    par(mar=c(0,5,2,1))
    plot(N[[2]], rep(2,length(N[[2]])), pch=124, col=col1, ylim=c(.5,2.5), xaxt='n', yaxt='n', ann=FALSE)
    points(N[[1]], rep(1,length(N[[1]])), pch=124, col=col2)
    mtext("Original spiketimes", las=1, side=3, adj=0, cex=.75)
    plot(s2, rep(2,length(s2)), pch=124, col=col1, ylim=c(.5,2.5), xaxt='n', yaxt='n', ann=FALSE)
    points(s1, rep(1,length(s1)), pch=124, col=col2)
    mtext("Rescaled spiketimes", las=1, side=3, adj=0, cex=.75)
dev.off()

plot(x,Lam[,1], type='l', col=col2)
lines(x,Lam[,2], type='l', col=col1)

layout(1:2)
hist(w1, breaks=seq(0,10,.1), prob=TRUE, col=col2, border=NA, ylim=c(0,1))
par(new=TRUE)
hist(w2, breaks=seq(0,10,.1), prob=TRUE, col=col1, border=NA, ylim=c(0,1))
curve(dexp, add=TRUE, from=1e-8, col=add.alpha('black',.75), lwd=2)

hist(v1, breaks=seq(0,10,.1), prob=TRUE, col=col2, border=NA, ylim=c(0,1))
par(new=TRUE)
hist(v2, breaks=seq(0,10,.1), prob=TRUE, col=col1, border=NA, ylim=c(0,1))
curve(dexp, add=TRUE, from=1e-8, col=add.alpha('black',.75), lwd=2)

par(mfrow=c(2,2))

ks.test(w1,'pexp')
ks.test(w2,'pexp')
ks.test(v1,'pexp')
ks.test(v2,'pexp')

pdf(paste0(plotlib, "embrechtsQQ1.pdf"), height=4, width=8)
    par(mfrow=c(1,2),mar=c(3,3,1,1), oma=c(0,0,0,0), bty='n')
    qqplot(w1, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col2, pch=16); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2)
    qqplot(w2, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col1, pch=16); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2)
dev.off()

pdf(paste0(plotlib, "embrechtsQQ2.pdf"), height=4, width=8)
    par(mfrow=c(1,2),mar=c(3,3,1,1), oma=c(0,0,0,0), bty='n')
    qqplot(v1, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col2, pch=16); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2)
    qqplot(v2, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col1, pch=16); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2)
dev.off()

layout(1)
pdf(paste0(plotlib, "embrechtsKS.pdf"), height=8, width=9)
    par(mfrow=c(2,2),mar=c(4,7,1,1), oma=c(0,0,0,0), bty='n')
    SSGLM::ks_plot(w1, col1 = col2, col2 = add.alpha('black',.15))
    mtext("Uncorrected rescaled waiting times w", 2, line=5.5)
    SSGLM::ks_plot(w2, col1 = col1, col2 = add.alpha('black',.15))
    
    SSGLM::ks_plot(v1, col1 = col2, col2 = add.alpha('black',.15))
    mtext("Mean corrected rescaled waiting times v", 2, line=5.5)
    SSGLM::ks_plot(v2, col1 = col1, col2 = add.alpha('black',.15))
dev.off()
    
poisson_lines <- function(s,ci=.95){
  a = qnorm((1+ci)/2)
  mu = length(s)/abs(diff(range(s)))
  mu.hi  = mu + a*sqrt(mu/length(s))
  mu.lo  = mu - a*sqrt(mu/length(s))
  abline(0,mu)
  abline((mu.hi-mu)*max(s),mu, lty=3)
  abline((mu.lo-mu)*max(s),mu, lty=3)
  return(c(mu,mu.lo,mu.hi))
}

# layout(1:2)
# plot(s1,1:length(s1), type='s');poisson_lines(s1);poisson_lines(s1,.99)
# plot(s2,1:length(s2), type='s');poisson_lines(s2);poisson_lines(s2,.99)


mu0  = n[1]
n0   = V[1,1]
B0   = d[1]
pp   = N[[1]]
Tmax = max(pp)
pars = c(mu0,n0,B0)

f = function(pars){
  -criticality::log_lik(pp, Tmax, pars[1], pars[2], pars[3])
}

system.time({
  res = constrOptim(c(mu0,n0,B0), f, ui = diag(3), ci=rep(0,3), method="Nelder-Mead")
})

pars
par.est = res$par


lam.est = criticality::cond_int(x, pp, par.est[1], par.est[2], par.est[3])
Lam.est = cumsum(lam.est)
L = data.frame(x=x,y=Lam.est)
s = interpol(as.matrix(L), x=pp)
w = diff(s)
v = w/mean(w);mean(w);sd(w)

pdf(paste0(plotlib,"estimated.pdf"), height=5, width=7)
    col3 = add.alpha('black',.5)
    par(mfrow=c(4,1), bty='n', oma=c(1,0,0,0), mar=c(2,5,1,1), las=1)
    layout(1:4,heights=c(8,8,2,2))
    tmp = as.Date(paste0(seq(1995,2010,5),"-01-01"))
    plot(pp,M[[1]], type='h', col=col3, ylim=c(-.1,0), xaxt='n', ann=FALSE)
    abline(h=seq(-.1,.1,.025), col=add.alpha('gray',.5))
    abline(h=0)
    axis(side = 1, at=as.numeric(tmp)-min(as.numeric(pfc$date[idx])), labels=substr(tmp,1,4))
    mtext("Excess log-returns", side=2, line=3.5, las=0, cex=.75)
    
    plot(x,lam[,1], type='l', col=col2, xaxt='n', ylim=c(0,.5), ann=FALSE)
    lines(x,lam.est, type='l', col=col3)
    axis(side = 1, at=as.numeric(tmp)-min(as.numeric(pfc$date[idx])), labels=substr(tmp,1,4))
    abline(h=seq(0,.5,.1), col=add.alpha('gray',.5))
    mtext("Estimated intensity", side=2, line=3.5, las=0, cex=.75)
    
    par(mar=c(0,5,2,1))
    plot(pp, rep(1,length(pp)), pch=124, col=col3, ylim=c(.5,2.5), xaxt='n', yaxt='n', ann=FALSE)
    mtext("Original spiketimes", las=1, side=3, adj=0, cex=.75)
    plot(s, rep(1,length(s)), pch=124, col=col3, ylim=c(.5,2.5), xaxt='n', yaxt='n', ann=FALSE)
    mtext("Rescaled spiketimes", las=1, side=3, adj=0, cex=.75)
dev.off()

pdf(paste0(plotlib,"diagnostics.pdf"), height=9, width=9)
    par(mfrow=c(2,2),mar=c(4,7,1,1), oma=c(0,0,3,0), bty='n')
    qqplot(w, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col3, pch=16, ann=FALSE); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2, las=0)
    mtext("Uncorrected rescaled waiting times w", 2, line=5.5, las=0)
    mtext("QQ plot", 3, line=2, las=0)
    SSGLM::ks_plot(w, col1 = col3, col2 = add.alpha('black',.15))
    mtext("KS plot", 3, line=2, las=0)
    
    qqplot(v, rexp(100, 1), xlim=c(0,8), ylim=c(0,8), col=col3, pch=16, ann=FALSE); abline(0,1)
    mtext("Theoretical", 1, 2)
    mtext("Observed", 2, 2, las=0)
    mtext("Mean corrected rescaled waiting times w", 2, line=5.5, las=0)
    SSGLM::ks_plot(v, col1 = col3, col2 = add.alpha('black',.15))
dev.off()

ks.test(w,'pexp')
ks.test(v,'pexp')
