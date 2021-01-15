library(bootha)
library(misc)
clean_up()
set.seed(1234)
layout(1)

pow_kernel <- function(x,a=.5,b=.75,k=.8){
  a*b/( (1+b*x)^(1+k) )
  
}
exp_kernel <- function(x,a=.5,b=.75){
  a*exp(-b*x)
}

a = .5
b = .75
k = .8

Tmax = 8000
pars = c(.10,.5,.75,.8)
pars[2]/pars[4]
system.time({
  pp = bootha::simulateHawkes(Tmax, pars)  
})


x  = seq(0,Tmax,1)
y  = intensity(x, pp, pars)
Nt = length(pp)

Nt
# layout(1:2)
# plot(x,y, type='l', ylim=c(-.1*max(y),max(y)))
# points(pp, rep(-.1*max(y),Nt), pch=124,col = add.alpha('black',.5))
# plot(x,integrated_intensity(x, pp, pars), type='l')


system.time({
  mle.pow = mle(pp, Tmax, pars, dt=1) # fit power kernel
  mle.exp = mle(pp, Tmax, pars[1:3], dt=0, ui=diag(c(1,-1,1)), ci=c(0,-1,0))   # fit exp kernel
})

L.pow = integrated_intensity(x, pp, mle.pow)
L.exp = integrated_intensity(x, pp, mle.exp)

s.pow = interpol(cbind(x,L.pow), pp)
s.exp = interpol(cbind(x,L.exp), pp)

layout(1)
plot(s.pow,s.exp, type='n')
abline(0,1, lty=3)
points(s.pow, s.exp, pch=16, col=add.alpha('red',.05))
mtext("Time transformed events: power law vs expontential kernel", side=3, line=2)


w.pow = diff(s.pow)
w.exp = diff(s.exp)

ks.test(w.pow,'pexp')
ks.test(w.exp,'pexp')


par(mfrow=c(3,2), bty='n')

plot(x,y, type='l', ylim=c(-.1*max(y),max(y)))
points(pp, rep(-.1*max(y),Nt), pch=124,col = add.alpha('black',.5))
mtext("Simulated power law process with intensity", side=3, line=2)

plot(x,integrated_intensity(x, pp, mle.pow), type='l')
lines(x,integrated_intensity(x, pp, mle.exp), col='red')
legend("topleft",c("Power law", "Exp kernel"), col=c('black','red'), lty=1, bty='n', lwd=3)
mtext("Integrated intensities", side=3, line=2)

plot(L.pow-L.exp, type='l')
mtext("Integrated intensity difference: powerlaw - exp kernel", side=3, line=2)

plot(s.pow-s.exp, type='l')
mtext("Time transformed events difference: powerlaw - exp kernel", side=3, line=2)

SSGLM::ks_plot(w.pow, col2=add.alpha('steelblue',.35))
mtext("Power law KS plot", side=3, line=2)
mtext(paste0("p-value ",round(ks.test(w.pow,'pexp')$p.value,3)), side=3, line=1, cex=.75)

SSGLM::ks_plot(w.exp, col2=add.alpha('steelblue',.35))
mtext("Exp kernel law KS plot", side=3, line=2)
mtext(paste0("p-value ",round(ks.test(w.exp,'pexp')$p.value,3)), side=3, line=1, cex=.75)

ks.test(w.pow,'pexp')
ks.test(w.exp,'pexp')

# fn = "~/Downloads/camilla.Rda"
# camilla = list(pp = pp, L.pow = L.pow, L.exp=L.exp, s.pow = s.pow, s.exp = s.exp)
# save(camilla, file=fn)
# str(camilla)


