library(ppboot)
library(misc)
misc::clean_up()

mu   = .75
a    = .5
b    = 1.1
Tmax = 100
dt   = 1e-2

# Test error difference of fb and rb (they should be equal...)
  set.seed(1234)
  hp   = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
  N    = length(hp)
  
  x = seq(0,Tmax,dt)
  l = data.frame(x=x,y=l.hawkes(x,hp, mu, a, b))
  L = data.frame(x=x,y=L.hawkes(x,hp, mu, a, b))
  
  s  = interpol(L,x=hp)
  fb = interpol(L,y=s)
  
  rb = NULL
  for(i in 1:length(s)){
    f   = L.hawkes(x,rb, mu, a, b)-s[i]
    tmp = interpol(cbind(x,f),y=0)
    rb  = c(rb,tmp)
    cat("\r",round(100*i/length(s),2),"% done!     ")
  }
  
  plot((rb-fb)^2, type='l')
  sum((rb-fb)^2)


# Visual inspection (use much shorter process!)
  Tmax = 10
  set.seed(1234)
  hp   = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
  N    = length(hp)
  
  x = seq(0,Tmax,dt)
  l = data.frame(x=x,y=l.hawkes(x,hp, mu, a, b))
  L = data.frame(x=x,y=L.hawkes(x,hp, mu, a, b))
  
  s  = interpol(L,x=hp)
  fb = interpol(L,y=s)
  
  rb = NULL
  for(i in 1:length(s)){
    f   = L.hawkes(x,rb, mu, a, b)-s[i]
    tmp = interpol(cbind(x,f),y=0)
    rb  = c(rb,tmp)
    cat("\r",round(100*i/length(s),2),"% done!     ")
  }
  
  Smax = max(L$y)

col1 = add.alpha('black',.95)
col2 = add.alpha('red',.75)
col3 = add.alpha('dodgerblue',.75)

# Check that fb and rb produce the same process with the same input Exp(1) variables and integrated conditional intensity
  layout(1:2)
  par(mar=c(3,3,1,1), bty='n', oma=c(0,0,0,0))
  plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,Smax), bty='n')
  lines(L$x, L$y, col=col1)
  lines(l$x, l$y, col=col1)
  
  points(hp, rep(0,N), pch=1, col=col1)
  points(rep(0,N), s, pch=16, col=col1)
  segments(hp,0,hp,s, lty=3, col=col1)
  segments(0,s,hp,s, lty=3, col=col1)
  
  
  points(fb, rep(0,N), pch=16, col=col2, cex=.75)
  segments(fb,0,fb,s, lty=3, col=col2)
  segments(0,s,fb,s, lty=3, col=col2)
  
  points(rb, rep(0,N), pch=16, col=col3, cex=.75)
  segments(rb,0,rb,s, lty=3, col=col3)
  segments(0,s,rb,s, lty=3, col=col3)



# Compare rb and fb with simulated Exp(1) variables

set.seed(12345)
s  = cumsum(rexp(2*N))
idx.fb = which(s<= max(L$y))
fb = interpol(L,y=s[idx.fb])

rb = NULL
for(i in 1:length(s)){
  f   = L.hawkes(x,rb, mu, a, b)-s[i]
  try({tmp = interpol(cbind(x,f),y=0)}, silent = TRUE)
  if(!is.na(tmp)){
    rb  = c(rb,tmp)  
  }
  cat("\r",round(100*i/length(s),2),"% done!     ")
}

l.rb = l.hawkes(x,rb, mu, a, b)
L.rb = L.hawkes(x,rb, mu, a, b)
l.fb = l.hawkes(x,fb, mu, a, b)
L.fb = L.hawkes(x,fb, mu, a, b)

plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,Smax), bty='n')
lines(L$x, L$y, col=add.alpha(col1,.4))
lines(L$x, l$y, col=add.alpha(col1,.4))

lines(x,l.fb, col=col2)
lines(x,L.fb, col=col2)
lines(x,l.rb, col=col3)
lines(x,L.rb, col=col3)


# Simulate a bunch of rb and fbs and compare
  layout(1)
  set.seed(12345)
  
  fb = rb = hps = list()
  B = 399
  for(j in 1:B){
    s  = cumsum(rexp(2*N))
    idx.fb = which(s<= max(L$y))
    fb[[j]] = interpol(L,y=s[idx.fb])
    rb.tmp = NULL
    for(i in 1:length(s)){
      f   = L.hawkes(x,rb.tmp, mu, a, b)-s[i]
      try({tmp = interpol(cbind(x,f),y=0)}, silent = TRUE)
      if(!is.na(tmp)){
        rb.tmp  = c(rb.tmp,tmp)  
      }
      
    }
    rb[[j]] = rb.tmp
    hps[[j]] = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
  
    cat("\r",round(100*j/B,2),"% done!     ")
  }


# Plot FIB and RB vs original
  plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,2*Tmax), bty='n')
  
  for(i in 1:B){
    l.fb = l.hawkes(x,fb[[i]], mu, a, b)
    L.fb = L.hawkes(x,fb[[i]], mu, a, b)
    l.rb = l.hawkes(x,rb[[i]], mu, a, b)
    L.rb = L.hawkes(x,rb[[i]], mu, a, b)
    # L.hp = L.hawkes(x,hps[[i]])
    
    # lines(x,l.fb, col=add.alpha(col2,.5))
    lines(x,L.fb, col=add.alpha(col2,.5))
    # lines(x,l.rb, col=add.alpha(col3,.5))
    lines(x,L.rb, col=add.alpha(col3,.5))
    
    # lines(x,L.rb, col=add.alpha(col1,.1))
  }
  
  lines(L$x, L$y, col=add.alpha(col1,.75), lwd=3)
  # lines(L$x, l$y, col=add.alpha(col1,.75), lwd=3)
  
  legend("topleft",c("Observed","FIB","RB"), lty=c(1,1,1), lwd=3, col=c(col1,col2,col3), bty='n')
  mtext("t",side=1, line=2)
  mtext(expression(Lambda~"(t)"),side=2, line=2)
  mtext(paste(B,"bootstrap samples"),side=3, line=0)

# Plot FIB and RB vs multiple simulated Hawkes (with similar pars)
  par(mfrow=c(1,3), mar=c(3,4,2,1), oma=c(0,0,0,0))
  col4 = add.alpha('orange',.75)
  plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,2*Tmax), bty='n')
  
  for(i in 1:B){
    # L.fb = L.hawkes(x,fb[[i]])
    # L.rb = L.hawkes(x,rb[[i]])
    L.hp = L.hawkes(x,hps[[i]])
    # lines(x,L.fb, col=add.alpha(col2,.5))
    # lines(x,L.rb, col=add.alpha(col3,.5))
    lines(x,L.hp, col=add.alpha(col1,.25))
  }
  lines(L$x, L$y, col=col4, lwd=3)

  
  legend("topleft",c("Simulated Hawkes","Observed Hawkes","FIB","RB"), lty=c(1,1,1), lwd=3, col=c(col1,col4,col2,col3), bty='n')
  mtext("t",side=1, line=2)
  mtext(expression(Lambda~"(t)"),side=2, line=2)
  

  plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,2*Tmax), bty='n')
  mtext(paste(B,"bootstrap samples"),side=3, line=0)
  
  for(i in 1:B){
    L.fb = L.hawkes(x,fb[[i]], mu, a, b)
    # L.rb = L.hawkes(x,rb[[i]])
    # L.hp = L.hawkes(x,hps[[i]])
    lines(x,L.fb, col=add.alpha(col2,.35))
    # lines(x,L.rb, col=add.alpha(col3,.5))
    # lines(x,L.hp, col=add.alpha(col1,.25))
  }
  lines(L$x, L$y, col=col4, lwd=3)
  
  plot(0,0, type='n', xlim=c(0,Tmax), ylim=c(0,2*Tmax), bty='n')
  
  for(i in 1:B){
    # L.fb = L.hawkes(x,fb[[i]])
    L.rb = L.hawkes(x,rb[[i]], mu, a, b)
    # L.hp = L.hawkes(x,hps[[i]])
    # lines(x,L.fb, col=add.alpha(col2,.5))
    lines(x,L.rb, col=add.alpha(col3,.35))
    # lines(x,L.hp, col=add.alpha(col1,.25))
  }
  lines(L$x, L$y, col=col4, lwd=3)
  
  
# Compare histograms of Ns (number of events) in each bootstrap process
  layout(1)

  hp.h = unlist(lapply(hps,length))
  rb.h = unlist(lapply(rb,length))
  fb.h = unlist(lapply(fb,length))
  
  hist(hp.h, breaks=seq(0,50,1), border=NA, col=add.alpha(col1,.3), prob=TRUE, ylim=c(0,.1))
  par(new=TRUE)
  hist(fb.h, breaks=seq(0,50,1), border=NA, col=add.alpha(col2,.3), prob=TRUE, ylim=c(0,.1))
  par(new=TRUE)
  hist(rb.h, breaks=seq(0,50,1), border=NA, col=add.alpha(col3,.3), prob=TRUE, ylim=c(0,.1))
  
  
  
  