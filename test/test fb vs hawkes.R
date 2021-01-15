library(ppboot)
library(misc)
misc::clean_up()

set.seed(1234)

Tmax = 10
mu   = .75
a    = .5
b    = 1.1
Tmax = 10
dt   = 1e-2
N    = 1000

hps = list()
for(i in 1:N){
  hps[[i]] = unlist(simulateHawkes(mu,a,b,Tmax))
}

x  = seq(0,Tmax,dt)
Ls = ls = matrix(NA, nr=length(x), nc=N)


for(i in 1:N){
  ok = FALSE
  try({ ltmp = make_l(hps[[i]],from = 0, to=Tmax, l.aux, dt = dt)
  ok = TRUE
  },silent=TRUE)
  if(ok){
    ls[,i] = ltmp$l
    Ls[,i] = ltmp$L
  }
}
unique(which(is.na(Ls),arr.ind = TRUE)[,2]) # Problematic columns...

i = 1000
tmp = hps[[i]]
ltmp = make_l(hps[[i]],from = 0, to=Tmax, l.aux, dt = dt)
ltmp
plot(ltmp$t, ltmp$l, type='l', ylim=c(0,10))
lines(ltmp$t, ltmp$L)
points(tmp, rep(0,length(tmp)))

y = interp.L(ltmp[,c(1,3)],x=tmp)
points(rep(0,length(y)),y)
hist(y, breaks=100)

idx = seq(1,nrow(Ls),1000)
# dim(Ls[idx,])
# layout(1:3)
# plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,max(Ls,na.rm = TRUE)), bty='n')
# for(i in 1:N){
#   lines(x[idx],Ls[idx,i], col=add.alpha('black',.2))
# }
# plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,max(ls,na.rm = TRUE)), bty='n')
# for(i in 1:N){
#   lines(x[idx],ls[idx,i], col=add.alpha('black',.2))
# }
Ns = unlist(lapply(hps, length))
# hist(Ns, breaks=seq(0,max(Ns),1), xlim=range(pretty(range(Ns))), ylim=c(0,.1), border=NA, col=add.alpha('black',.4), prob=TRUE)


S = rb = list()
for(i in 1:100){
  S[[i]]  = cumsum(rexp(1.5*max(Ns)))
  try({
    # int     = make_l(S[[i]],0,Tmax,l.aux,1e-3)
    # fb[[i]] = invert.L(S[[i]],int[,c(1,3)],Tmax, dt=1e-3) 
    rb[[i]] = invert.L(S[[i]],L.aux,Tmax, dt=1e-3)   
  })
  cat("\r",i)
}

# Ss = unlist(lapply(rb,nrow))
# par(new=TRUE)
# hist(Ss, breaks=seq(0,max(Ss),1), xlim=range(pretty(range(Ns))), ylim=c(0,.05), border=NA, col=add.alpha('red',.4), prob=TRUE)
# range(Ss)
# range(Ns)
layout(1)
# plot(x,ls[,1], type='l')

i=1

rb.Ls = rb.ls = matrix(NA, nr=length(x), nc=length(rb))
for(i in 1:length(rb)){
  ok = FALSE
  rb.tmp = rb[[i]]
  try({ ltmp = make_l(rb.tmp$t,from = 0, to=Tmax, l.aux, dt = dt)
  ok = TRUE
  },silent=TRUE)
  if(ok){
    rb.ls[,i] = ltmp$l
    rb.Ls[,i] = ltmp$L
  }
}
unique(which(is.na(rb.Ls),arr.ind = TRUE)[,2]) # Problematic columns...

idx = seq(1,nrow(Ls),1000)
layout(1:2)
plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,max(Ls,na.rm = TRUE)), bty='n')
for(i in 1:N){
  lines(x[idx],Ls[idx,i], col=add.alpha('black',.2))
}
for(i in 1:length(rb)){
  lines(x[idx],rb.Ls[idx,i], col=add.alpha('red',.2))
}

plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,max(ls,na.rm = TRUE)), bty='n')
for(i in 1:N){
  lines(x[idx],ls[idx,i], col=add.alpha('black',.2))
}
for(i in 1:length(rb)){
  lines(x[idx],rb.ls[idx,i], col=add.alpha('red',.2))
}

