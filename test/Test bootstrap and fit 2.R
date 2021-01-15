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

library(emhawkes)

mu.in = length(hp)/Tmax   # Set baseline guess to observed events/interval
a.in = mu.in/10           # Guess intensity jump is 10% of baseline
b.in = 2*a.in             # Ensure a/b ratio is less than 1

hspec <- new("hspec", mu=mu.in, alpha=a.in, beta=b.in)
invisible(capture.output({mle = hfit(hspec, inter_arrival = diff(hp), lambda0 = mu.in)}))  

mu.in = mle$estimate[1]
a.in = mle$estimate[2]
b.in = mle$estimate[3]
bootspec <- new("hspec", mu=mu.in, alpha=a.in, beta=b.in)

l = data.frame(x=x,y=l_exp_kernel(x,hp,mu.in,a.in,b.in))
L = data.frame(x=x,y=L_exp_kernel(x,hp,mu.in,a.in,b.in))

s  = interpol(as.matrix(L),x=hp)

kern <- function(x,y){
  if(is.null(y)){
    y = numeric(0)
  }
  out = L_exp_kernel(x,y,mu.in,a.in,b.in)
}


B    = 399
Smax = max(L$y)

system.time({
  fb = pp.boot(L,Smax,kern,B,NULL,parallel = TRUE)
})


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
col1 = add.alpha('dodgerblue',.5)
col2 = add.alpha('dodgerblue',.5)
col3 = add.alpha('dodgerblue',.5)

mle.est = mle$estimate
mle.se = sqrt(diag(vcov(mle)))


range(mus)
hist(mus, breaks=seq(-20,20,.01), xlim=c(0,2), border=NA, col=col1, prob=TRUE); abline(v=mu,lty=3)
ses = (out2[1,3,]-out2[1,2,])/1.96
sds = mean(ses, na.rm = TRUE)+seq(-sd(ses, na.rm = TRUE),sd(ses, na.rm = TRUE),length=10)
mns = mean(mus, na.rm = TRUE)
# for(i in 1:10){
#   curve(dnorm(x,mns,sds[i]), add=TRUE, lty=3, col=add.alpha('black',.5))
# }
curve(dnorm(x,mns,sd(mus)), add=TRUE, lty=1, col=add.alpha('red',.75))
curve(dnorm(x,mle.est[1], mle.se[1]), add=TRUE, lty=1, col=add.alpha('orange',.75))

plot(mus, col=col1, pch=16, cex=.75, ylim=c(0,2)); abline(h=mu,lty=3)


range(as)
hist(as, breaks=seq(-20,20,.01), xlim=c(0,2), border=NA, col=col2, prob=TRUE); abline(v=a,lty=3)
ses = (out2[2,3,]-out2[2,2,])/1.96
sds = mean(ses, na.rm = TRUE)+seq(-sd(ses, na.rm = TRUE),sd(ses, na.rm = TRUE),length=10)
mns = mean(as, na.rm = TRUE)
# for(i in 1:10){
#   curve(dnorm(x,mns,sds[i]), add=TRUE, lty=3, col=add.alpha('black',.5))
# }
curve(dnorm(x,mns,sd(as)), add=TRUE, lty=1, col=add.alpha('red',.75))
curve(dnorm(x,mle.est[2], mle.se[2]), add=TRUE, lty=1, col=add.alpha('orange',.75))

plot(as, col=col2, pch=16, cex=.75, ylim=c(0,2)); abline(h=a,lty=3)

range(bs)
hist(bs, breaks=seq(-20,100,.01), xlim=c(0,2), border=NA, col=col3, prob=TRUE); abline(v=b,lty=3)
ses = (out2[3,3,]-out2[3,2,])/1.96

sds = mean(ses, na.rm = TRUE)+seq(-sd(ses[ses<10], na.rm = TRUE),sd(ses[ses<10], na.rm = TRUE),length=10)
mns = median(bs, na.rm = TRUE)
# for(i in 1:10){
#   if(sds[i]>0){
#     curve(dnorm(x,mns,sds[i]), add=TRUE, lty=3, col=add.alpha('black',.5))  
#   }
# }
curve(dnorm(x,mns,sd(bs[bs<5])), add=TRUE, lty=1, col=add.alpha('red',.75))
curve(dnorm(x,mle.est[3], mle.se[3]), add=TRUE, lty=1, col=add.alpha('orange',.75))
plot(bs, col=col3, pch=16, cex=.75, ylim=c(0,2)); abline(h=b,lty=3)


layout(1:2)
hist(as/bs, col=col1, border=NA, xlim=c(0,1), breaks=seq(-10,20,.01), prob=TRUE)
abline(v=a/b, lty=3)
plot(as/bs, col=col1, pch=16, cex=.75, ylim=c(0,1))
abline(h=a/b, lty=3)

filen = "/Users/jmf408/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/rb_test_par.Rda"
save(out,file=filen)





# rgl::plot3d(mus,as,bs, xlim=c(0,2), ylim=c(0,2), zlim=c(0,2), col='dodgerblue')
# rgl::points3d(mu,a,b, col='red', size=10)
# rgl::points3d(mle.est[1],mle.est[2],mle.est[3], col='orange', size=10)
# round(mle.est,3)
# round(1.96*mle.se,3)

