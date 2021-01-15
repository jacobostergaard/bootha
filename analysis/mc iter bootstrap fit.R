library(emhawkes)
library(hawkes)
library(ppboot)
library(misc)
library(parallel)
clean_up()

set.seed(1234)



n_iter  = 500

B    = 299
mu   = 1
a    = 4
b    = 5
Tmax = 20

parametric = TRUE
boot.type = 'fib'


boot <- function(obs, B, pars, parametric = TRUE, type='fib', dt=1e-2){
  
  mu.in = pars[1]
  a.in = pars[2]
  b.in = pars[3]
  
  x = seq(0,Tmax,dt)
  l = data.frame(x=x,y=l_exp_kernel(x,obs,mu.in,a.in,b.in))
  L = data.frame(x=x,y=L_exp_kernel(x,obs,mu.in,a.in,b.in))
  s = interpol(as.matrix(L),x=obs)
  if(type=='fib'){
    kern = NULL
  } else{
    kern <- function(x,y){
      if(is.null(y)){
        y = numeric(0)
      }
      out = L_exp_kernel(x,y,mu.in,a.in,b.in)
    }  
  }
  if(parametric){
    s.events = NULL
  } else{
    s.events = s
  }
  
  Smax = max(L$y)
  bb = pp.boot(L = L, Smax = Smax, kernel = kern, B = B, s.events = s.events, parallel = TRUE)
  
  return(bb) 
}
fit_boot <- function(bb, pars){
  spec = new("hspec", mu=pars[1], alpha=pars[2], beta=pars[3])
  # tmp.fit = lapply(bb, function(x) fit = quiet(emhawkes::hfit(spec, diff(c(0,x)), lambda0=pars[1]))
  tmp.fit = parallel::mclapply(bb, function(x) fit = quiet(emhawkes::hfit(spec, diff(c(0,x)), lambda0=pars[1])), mc.cores = 4)
  tmp.est = matrix(unlist(lapply(tmp.fit, function(x) x$estimate)), nc=3, byrow=TRUE)
  tmp.est = cbind(tmp.est,unlist(lapply(tmp.fit, function(x) x$maximum)))
  colnames(tmp.est) = c("mu","alpha","beta","loglik")
  rownames(tmp.est) = 1:length(bb)
  return(as.data.frame(tmp.est))
}
mc <- function(n_iter, mu, a, b, Tmax, B, parametric = TRUE, boot.type = 'fib', verbose=TRUE, filename=NULL){
  
  tmr = numeric(0) # Timer
  boot.pval = chi2.pval = numeric(n_iter)
  
  for(i in 1:n_iter){
    tic = Sys.time()
    # obs = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
    spec = new("hspec", mu=mu, alpha=a, beta=b)
    obs = hsim(spec, size=10*mu/(1-a/b)*Tmax, lambda0 = mu)$arrival[-1]
    obs = obs[obs<=Tmax]
    mle = quiet(emhawkes::hfit(spec, diff(c(0,obs))))
    
    
    theta0 = c(mu,a,b)
    thetah = mle$estimate
    
    lrt = 2*(llk(thetah, obs)-llk(theta0, obs))
    
    # print("Restricted parametric FIB")
    # system.time({ 
      boot.res = boot(obs,B=B, theta0, parametric = parametric, type = boot.type) 
      # })
    # print("Fitting bootstrap MLEs")
    # system.time({ 
      boot.fit = fit_boot(boot.res, theta0) 
      # })
    
    boot.llk0 = unlist(lapply(boot.res, function(x) llk(theta = theta0, x)))
    boot.llkh = boot.fit$loglik
    boot.lrt  = 2*(boot.llkh - boot.llk0)
    boot.pval[i] = mean(lrt<boot.lrt)
    chi2.pval[i] = 1-pchisq(lrt,3)
    
    # Timer part, displays approx. time left
    toc = Sys.time()
    tmr = c(tmr,difftime(toc,tic, units = "secs"))
    t_left = mean(tmr)*(n_iter-i)
    pct = round(100*i/n_iter)
    if(verbose){
      cat("\rProgress:",pct,"%  Time left approx.", t_left,"seconds")  
    }
  
    if(!is.null(filename)){
      boot = list(boot = boot.pval[1:i], chi2 = chi2.pval[1:i])  
      save(boot, file = filename)
    }
    
  }  
  
  out = list(boot = boot.pval, chi2 = chi2.pval)
  
  return(out)
}


tic = Sys.time()
fn  = "/Users/jmf408/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/ppboot/data/bootres.Rda"
res = mc(n_iter, mu, a, b, Tmax, B, parametric, boot.type, filename=fn)

toc = Sys.time()
cat("\nTotal time:",as.numeric(difftime(toc,tic, units = "secs")))

# load(fn)
# lapply(boot, function(x) mean(x<0.05))
