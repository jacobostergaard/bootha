library(emhawkes)
library(hawkes)
library(ppboot)
library(misc)
library(parallel)
clean_up()


mu    = 1
a    = 4
b    = 5
Tmax = 100

spec <- new("hspec", mu=mu, alpha=a, beta=b)

set.seed(1234)

MCiter = 1000

f <- function(i=1){
  obs = c(0,unlist(hawkes::simulateHawkes(mu,a,b,Tmax)))
  fit = quiet(emhawkes::hfit(spec, diff(obs)))
  est = fit$estimate
  ll  = -hawkes::likelihoodHawkes(est[1],est[2],est[3],obs)  
  
  out = list(obs=obs,fit=fit,ll=ll)
  
  return(out)
}

fdir = "~/Library/Mobile Documents/com~apple~CloudDocs/UCPH/CIBS/Hawkes project/data/mc/"
# system.time({res = mclapply(1:10, f, mc.cores = 4)})

iter = 100
mc = 10
tmr = numeric(0) # Timer
for(i in 1:iter){
  tic = Sys.time()
  res = quiet(mclapply(1:mc, f, mc.cores = 3))
  sink()
  fn = paste0(fdir,"mc",i,".rda")
  save(res,file = fn)
  
  toc = Sys.time()
  tmr = c(tmr,difftime(toc,tic, units = "secs"))
  t_left = mean(tmr)*(iter-i)
  
  pct = round(100*i/iter)
  cat("\rProgress:",pct,"%  Time left approx.", t_left,"seconds")
}

allres = list()
k = 1
for(i in 1:iter){
  fn = paste0(fdir,"mc",i,".rda")
  load(fn)
  for(j in 1:mc){
    allres[[k]] = res[[j]]  
    k = k+1
  }
  
  cat("\ri=",i)
}

unlist(lapply(allres,function(x) length(x)))

res = allres[unlist(lapply(allres,length))>1]
obs = lapply(res, function(x) x$obs)
ll0 = unlist(lapply(obs, function(x) -hawkes::likelihoodHawkes(mu,a,b,x)))
lls = unlist(lapply(res, function(x) x$ll))

lrt = -2*(ll0-lls)
hist(lrt, breaks=100, prob=TRUE, border=NA, col=add.alpha('red',.45))
curve(dchisq(x,3), add=TRUE, lwd=2)


which(unlist(lapply(allres,length))==1) 
res = allres[]

abline(v = quantile(lrt,.95), col='red')
abline(v = qchisq(.95,3), lty=3)
