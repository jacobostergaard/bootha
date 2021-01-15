library(misc)
library(bootha)
clean_up()
set.seed(1234)

Tmax = 6000
pars = c(0.018, 0.74, 0.021) # Embrechts parameters

pp = simulateHawkes(Tmax, pars)

length(pp) # simulated observations

n0  = pars[2]
mu0 = n0*length(pp)/Tmax
B0  = pars[3]

# These constraints ensures mu>0, n<1 and B>0:
# ui0 = diag(c(1,-1,1))
# ci0 = c(0,-1,0)
(res = mle(pp, Tmax, pars))

log_lik(pp, Tmax, pars)
log_lik(pp, Tmax, res)





# boot = list()
# system.time({
#   boot[[1]] = bootstrap_hawkes1(pp, pars, Tmax, "fib", parametric=TRUE)
#   boot[[2]] = bootstrap_hawkes1(pp, pars, Tmax, "fib", parametric=FALSE)
#   boot[[3]] = bootstrap_hawkes1(pp, pars, Tmax, "rib", parametric=TRUE)
#   boot[[4]] = bootstrap_hawkes1(pp, pars, Tmax, "rib", parametric=FALSE)
# })

# These constraints ensures mu>0, n<1 and B>0:
ui0 = diag(c(1,-1,1))
ci0 = c(0,-1,0)


# These constraints ensures mu>0, n>0 and B>0:
# ui0 = diag(c(1,1,1))
# ci0 = c(0,0,0)
# for(b in 1:length(boot)){
#   print(mle(boot[[b]], Tmax, pars))
# }
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
fn = paste0(datalib,"boot_test.Rda")

B = 399
boot = list()
f <- function(b){
  pp = bootstrap_hawkes1(pp, pars, Tmax, "fib", parametric=TRUE)
  res = mle(pp, Tmax, pars)
  return(res)
}
res = matrix(unlist(lapply(1:B, f)),nc=3, byrow=TRUE)
colnames(res) = c("mu","n","B")
boot[[1]] = res
save(boot,file=fn)
notify("Parameteric fib done!")

f <- function(b){
  pp = bootstrap_hawkes1(pp, pars, Tmax, "fib", parametric=FALSE)
  res = mle(pp, Tmax, pars)
  return(res)
}
res = matrix(unlist(lapply(1:B, f)),nc=3, byrow=TRUE)
colnames(res) = c("mu","n","B")
boot[[2]] = res
save(boot,file=fn)
notify("Non-parameteric fib done!")

f <- function(b){
  pp = bootstrap_hawkes1(pp, pars, Tmax, "rib", parametric=TRUE)
  res = mle(pp, Tmax, pars)
  return(res)
}
res = matrix(unlist(lapply(1:B, f)),nc=3, byrow=TRUE)
colnames(res) = c("mu","n","B")
boot[[3]] = res
save(boot,file=fn)
notify("Parameteric rib done!")

f <- function(b){
  pp = bootstrap_hawkes1(pp, pars, Tmax, "rib", parametric=FALSE)
  res = mle(pp, Tmax, pars)
  return(res)
}
res = matrix(unlist(lapply(1:B, f)),nc=3, byrow=TRUE)
colnames(res) = c("mu","n","B")
boot[[4]] = res
save(boot,file=fn)
notify("Non-parameteric rib done!")

# B = 399
# library(foreach)
# 
# result <- foreach(x = 1:B) %do% f(x)
# parallel::mclapply(1:2, f, mc.cores=2)
# tmp.fit = parallel::mclapply(bb, function(x) fit = quiet(emhawkes::hfit(spec, diff(c(0,x)), lambda0=pars[1])), mc.cores = 4)
  
names(boot) = c("par_fib","npar_fib","par_rib","npar_rib")

boot_par_ci <- function(boot_in){
  tmp = t(rbind(apply(boot_in,2,mean),apply(boot_in,2,mean)-1.96*apply(boot_in,2,sd),apply(boot_in,2,mean)+1.96*apply(boot_in,2,sd)))
  colnames(tmp) = c("est","ci.lo","ci.hi")
  return(tmp)
}

pars
boot_par_ci(boot$par_fib)
boot_par_ci(boot$npar_fib)
boot_par_ci(boot$par_rib)
boot_par_ci(boot$npar_rib)

