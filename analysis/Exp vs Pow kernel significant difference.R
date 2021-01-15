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

sim_compare <- function(pars, Tmax){
  
    # system.time({
      pp = bootha::simulateHawkes(Tmax, pars)  
    # })
    
    x  = seq(0,Tmax,1)
    y  = intensity(x, pp, pars)
    Nt = length(pp)
    cat("\n",Nt," observations",sep="")
    # system.time({
      mle.pow = mle(pp, Tmax, pars, dt=1) # fit power kernel
      mle.exp = mle(pp, Tmax, pars[1:3], dt=0, ui=diag(c(1,-1,1)), ci=c(0,-1,0))   # fit exp kernel
    # })
    
    L.pow = integrated_intensity(x, pp, mle.pow)
    L.exp = integrated_intensity(x, pp, mle.exp)
    s.pow = interpol(cbind(x,L.pow), pp)
    s.exp = interpol(cbind(x,L.exp), pp)
    w.pow = diff(s.pow)
    w.exp = diff(s.exp)
    
    w = data.frame(pow=w.pow, exp=w.exp)
    print(ks.test(w$pow,'pexp'))
    print(ks.test(w$exp,'pexp'))

    layout(1:2)
    ks_plot(w$pow)
    ks_plot(w$exp)
    
    return(w)
}


Tmax = 4000

cat("\n---------------------------------------------------------------------------------------")

pars = c(.10,.5,.75,.8)
# pars[2]/pars[4]

cat("\nmu=",pars[1], sep="")
system.time({
  w = sim_compare(pars, Tmax)  
})
notify("Simulation done!")

cat("\n---------------------------------------------------------------------------------------")

pars = c(.20,.5,.75,.8)
# pars[2]/pars[4]

cat("\nmu=",pars[1], sep="")
system.time({
  w = sim_compare(pars, Tmax)  
})
notify("Simulation done!")


cat("\n---------------------------------------------------------------------------------------")

pars = c(.40,.5,.75,.8)
# pars[2]/pars[4]

cat("\nmu=",pars[1], sep="")
system.time({
  w = sim_compare(pars, Tmax)  
})
notify("Simulation done!")


cat("\n---------------------------------------------------------------------------------------")

pars = c(.80,.5,.75,.8)
# pars[2]/pars[4]

cat("\nmu=",pars[1], sep="")
system.time({
  w = sim_compare(pars, Tmax)  
})
notify("Simulation done!")


