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

mle  = mle(pp, Tmax, pars, dt=0, grad = function(x) exp_gradient(pp, Tmax, x), hessian=FALSE)


cat(paste0("\nBegin at: ", Sys.time()),"\n")
fn = paste0(datalib,"dji_boot.Rda")
B = 399
ncores = 1
dji.boot = list() 
system.time({
  dji.boot[[1]] = bootstrap_hawkes(pp, mle, Tmax, B, type='fib', parametric=TRUE, ncores=ncores)
})
system.time({
  dji.boot[[2]] = bootstrap_hawkes(pp, mle, Tmax, B, type='fib', parametric=FALSE, ncores=ncores)
})
system.time({
  dji.boot[[3]] = bootstrap_hawkes(pp, mle, Tmax, B, type='rib', parametric=TRUE, ncores=ncores)
})
system.time({
  dji.boot[[4]] = bootstrap_hawkes(pp, mle, Tmax, B, type='rib', parametric=FALSE, ncores=ncores)
})
names(dji.boot) = c("fib_par","fib_npar", "rib_par","rib_npar")

cat(paste0("\nDone at: ", Sys.time()))

# msg = paste0("DJI bootstrap done")
# notify(msg)
# cat("\r",msg, sep="")

# par(mfrow=c(2,2))
# for(i in 1:4){
#   plot(0,0, type='n', xlim=c(0,6500), ylim=c(0,600))
#   invisible(lapply(dji.boot[[i]], function(x) lines(x, 1:length(x), type='s', col=add.alpha('black',.3))))
# }

save(dji.boot, file=fn)

