library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

plotlib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/plots/"
datalib = "~/Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/bootha/data/"
load(paste0(datalib,"boot_embrechts5.Rda"))
load(paste0(datalib,"dowjones.Rda"))

pp   = dowjones$pp.neg
pars = dowjones$pars
Tmax = dowjones$Tmax

# These constraints ensures mu>0, n<1 and B>0:
# ui0 = diag(c(1,-1,1))
# ci0 = c(0,-1,0)
(res = mle(pp, Tmax, pars))
# 
# log_lik(pp, Tmax, pars)
# log_lik(pp, Tmax, res)

x = seq(0,Tmax,1)
L.exp = integrated_intensity(x, pp, res)
s.exp = interpol(cbind(x,L.exp), pp)

s = s.exp
w = diff(s)
v = w/mean(w);mean(w);sd(w)

layout(1)
acf(diff(pp))

par(mfrow=c(2,2))

plot_acf(w, lag=100)
plot_acf(v, lag=100)
plot_acf(w^2, lag=100)
plot_acf(v^2, lag=100)



# These constraints ensures mu>0, n<1 and B>0:
ui0 = diag(c(1,-1,1))
ci0 = c(0,-1,0)

# These constraints ensures mu>0, n>0 and B>0:
# ui0 = diag(c(1,1,1))
# ci0 = c(0,0,0)
# for(b in 1:length(boot)){
#   print(mle(boot[[b]], Tmax, pars))
# }
fn = paste0(datalib,"boot_embrechts.Rda")

booter <- function(B, type, parametric, dt=1, ncores=6){
  tmp = bootstrap_hawkes(pp, pars, Tmax, B, type, parametric,  dt, ncores)
  cat("    Bootstrap done, estimating mle...")
  res = parallel::mclapply(tmp, function(pp) mle(pp, Tmax, pars, dt=0), mc.cores=ncores)
  # res = lapply(tmp, function(pp) mle(pp, Tmax, pars) )
  res = matrix(unlist(res),nc=3, byrow=TRUE)
  colnames(res) = c("mu","n","B")
  return(res)
}

# Estimate time to completion
# B = 1
# system.time({
#   bootstrap_hawkes(pp, pars, Tmax, B, type='rib', parametric=TRUE,  dt=10, ncores=6)
# })
# 
# system.time({
#   bootstrap_hawkes(pp, pars, Tmax, B, type='rib', parametric=TRUE,  dt=5, ncores=6)
# })
# 
# system.time({
#   bootstrap_hawkes(pp, pars, Tmax, B, type='rib', parametric=TRUE,  dt=1, ncores=6)
# })

# dt = 1  => 399 rib samples in about  5h
# dt = 5  => 399 rib samples in about  1h
# dt = 10 => 399 rib samples in about .5h

fn = paste0(datalib,"boot_embrechts.Rda")

B = 399
boot = list()
cat(paste0("\nBegin at: ", Sys.time()),"\n")
boot[[1]] = booter(B, 'fib', parametric=TRUE); save(boot,file=fn); notify("Parameteric fib done!")
cat(paste0("\nPar fib done at: ", Sys.time()),"\n")
boot[[2]] = booter(B, 'fib', parametric=FALSE); save(boot,file=fn); notify("Non-parameteric fib done!")
cat(paste0("\nNon-par fib done at: ", Sys.time()),"\n")
boot[[3]] = booter(B, 'rib', parametric=TRUE, dt=1); save(boot,file=fn); notify("Parameteric rib done!")
cat(paste0("\nPar rib done at: ", Sys.time()),"\n")
boot[[4]] = booter(B, 'rib', parametric=FALSE, dt=1); save(boot,file=fn); notify("Non-parameteric rib done!")
cat(paste0("\nNon-par fib done at: ", Sys.time()))

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



