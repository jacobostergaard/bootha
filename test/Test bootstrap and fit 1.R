library(ppboot)
library(misc)
misc::clean_up()

mu   = .75
a    = .5
b    = 1.1
Tmax = 100
dt   = 1e-1

# Test error difference of fb and rb (they should be equal...)
set.seed(1234)
hp   = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
N    = length(hp)

x = seq(0,Tmax,dt)

l = data.frame(x=x,y=l_exp_kernel(x,hp,mu,a,b))
L = data.frame(x=x,y=L_exp_kernel(x,hp,mu,a,b))


s  = interpol(as.matrix(L),x=hp)

kern <- function(x,y){
  if(is.null(y)){
    y = numeric(0)
  }
  out = L_exp_kernel(x,y,mu,a,b)
}


B    = 399
Smax = max(L$y)

system.time({
  fb = pp.boot(L,Smax,NULL,B,s,parallel = TRUE)
})

library(emhawkes)

mu.in = length(hp)/Tmax   # Set baseline guess to observed events/interval
a.in = mu.in/10           # Guess intensity jump is 10% of baseline
b.in = 2*a.in             # Ensure a/b ratio is less than 1

hspec <- new("hspec", mu=mu.in, alpha=a.in, beta=b.in)
invisible(capture.output({mle = hfit(hspec, inter_arrival = diff(hp), lambda0 = mu.in)}))  

bootspec <- new("hspec", mu=mle$estimate[1], alpha=mle$estimate[2], beta=mle$estimate[3])
B = 5
res = array(dim=c(3,3,B), dimnames = list(c("mu","alpha","beta"), c("lo","est","hi")))
tmr = numeric(0) # Timer
system.time({
  for(i in 1:B){
    tic = Sys.time()
    invisible(capture.output({mle.tmp = hfit(bootspec, inter_arrival = diff(fb[[i]]), lambda0 = mle$estimate[1])}))  
    est = mle.tmp$estimate
    se  = sqrt(diag(vcov(mle.tmp)))
    tmp = cbind(est-1.96*se,est,est+1.96*se)
    res[,,i] = tmp
  
    toc = Sys.time()
    tmr = c(tmr,difftime(toc,tic, units = "secs"))
    t_left = mean(tmr)*(B-i)
    pct = round(100*i/B)
    cat("\rProgress:",pct,"%  Time left approx.", t_left,"seconds")
  } 
  cat("\n")
})

f <- function(i=1){
  tic = Sys.time()
  invisible(capture.output({mle.tmp = hfit(bootspec, inter_arrival = diff(fb[[i]]), lambda0 = mle$estimate[1])}))  
  est = mle.tmp$estimate
  se  = sqrt(diag(vcov(mle.tmp)))
  tmp = data.frame(lo = est-1.96*se,est = est, hi= est+1.96*se)
  return(tmp)
}

system.time({out = mclapply(1:B,f,mc.cores = 4)})

B = 399
system.time({out = mclapply(1:B,f,mc.cores = 4)})


notify("Boostrap done!")


out2 = array(unlist(out),dim=c(3,3,B), dimnames = list(c("mu","alpha","beta"), c("lo","est","hi")) )
mus = out2[1,2,]
as = out2[2,2,]
bs = out2[3,2,]

par(mfrow=c(3,2), bty='n')
col1 = add.alpha('red',.5)
col2 = add.alpha('red',.5)
col3 = add.alpha('red',.5)
hist(mus, breaks=seq(-20,20,.01), xlim=c(0,2), border=NA, col=col1, prob=TRUE); abline(v=mu,lty=3)
plot(mus, col=col1, pch=16, cex=.75, ylim=c(0,2)); abline(h=mu,lty=3)

hist(as, breaks=seq(-20,20,.01), xlim=c(-1,1), border=NA, col=col2, prob=TRUE); abline(v=a,lty=3)
plot(as, col=col2, pch=16, cex=.75, ylim=c(-1,1)); abline(h=a,lty=3)

hist(bs, breaks=seq(-20,20,.01), xlim=c(0,2), border=NA, col=col3, prob=TRUE); abline(v=b,lty=3)
plot(bs, col=col3, pch=16, cex=.75, ylim=c(0,2)); abline(h=b,lty=3)

layout(1)
plot(as/bs, col=col1, pch=16, cex=.75, ylim=c(-2,2))
abline(h=a/b, lty=3)

filen = "/Users/jmf408/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/fb_test.Rda"
save(out,file=filen)
