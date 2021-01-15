library(bootha)
library(misc)
misc::clean_up()
set.seed(1234)

hawkes_gradient <- function(Tmax, Nt, pars){
  u = pars[1]
  y = pars[2]
  b = pars[3]
  
  n = length(Nt)
  A = numeric(n)
  B = numeric(n)
  C = numeric(n)
  
  for(i in 2:n){
    ti    = Nt[i]-Nt[i-1]
    A[i] = exp( -b*ti )*( 1+A[i-1] )
    B[i] = ti*A[i]+exp(-b*ti)*B[i-1]
    C[i] = ti*(B[i]+exp(-b*ti)*B[i-1])+exp(-b*ti)*C[i-1]
    
  }
  
  Ttmp = Tmax-Nt
  dldu = -Tmax+sum( 1/(u+y*b*A) )
  dldy = sum(  exp(-b*Ttmp)-1 + b*A/(u+y*b*A)  )
  dldb = -y*sum( exp(-b*Ttmp)*Ttmp  ) + sum( (y*A-y*b*B)/(u+y*b*A) )
  
  out = c(dldu, dldy, dldb)
  return(out)
}
hawkes_hessian <- function(pp, Tmax, numerical=TRUE){
  # Find Hessian by numerical solution:
  f = function(pars){
    -log_lik(pp, Tmax, pars, 0)
  }
  g = function(pars){
    hawkes_gradient(Tmax, pp, pars)  
  }
  ui = diag(length(pars))
  ci = rep(0,length(pars))
  # If gradient can be supplied, it is possible to obtain the Hessian
  res = constrOptim(pars, f=f, grad = g, ui = ui, ci=ci, method="Nelder-Mead", hessian=TRUE)
  # names(mle) = c("mu","n","B")
  mle = res$par
  
  if(numerical){
    H   = res$hessian
  } else{
    
    A = ABC(mle,pp)$A
    B = ABC(mle,pp)$B
    C = ABC(mle,pp)$C
    
    tmp = Tmax-pp
    u   = mle[1]
    y   = mle[2]
    b   = mle[3]
    Huu = -sum((u+y*b*A)^-2)
    Hyy = -sum( (b*A/(u+y*b*A))^2 )
    Hbb = y*sum( exp(-b*tmp)*tmp^2 )+ sum( (y*b*C-2*y*B)/(u+y*b*A) ) - sum( (y*A-y*b*B)^2/(u+y*b*A)^2 )
    Huy = -sum( b*A/(u+y*b*A)^2 )

    Hub = -sum( (y*A-y*b*B)/(u+y*b*A)^2 )
    Hyb = -sum( tmp*exp(-b*tmp) ) + sum( (A-b*B)/(u+y*b*A) ) - sum( y*b*A*(A-b*B)/(u+y*b*A)^2 )
    H   = matrix(c(Huu,Huy,Hub,Huy,Hyy,Hyb,Hub,Hyb,Hbb), nr=3, nc=3, byrow=TRUE)
  }
  out = list(I=-solve(H), H=H, mle=mle)
  
  return(out) 
}   
plot_gradient <- function(par_no,pars,pp,Tmax, lo=.5, hi=1.5){
  # par_no = 3
  # pars = res
  # pp
  # Tmax
  # lo
  # hi
  
  nx   = 1000
  xmin = lo*pars[par_no]
  xmax = hi*pars[par_no]
  x = seq(xmin,xmax,length=nx)
  z = numeric(nx)
  for(i in 1:nx){
    if(par_no==1){
      z[i] = log_lik(pp,Tmax, c(x[i],pars[2],pars[3]))-log_lik(pp, Tmax, res)  
    } else if(par_no==2){
      z[i] = log_lik(pp,Tmax, c(pars[1],x[i],pars[3]))-log_lik(pp, Tmax, res)  
    } else if(par_no==3){
      z[i] = log_lik(pp,Tmax, c(pars[1],pars[2],x[i]))-log_lik(pp, Tmax, res)  
    }
  }
  plot(x,z, xaxt='n', type='l', col=add.alpha('black',.75), yaxt='n', bty='n', ann=FALSE)
  abline(-hgrad[par_no]*pars[par_no],hgrad[par_no], col=add.alpha('black',.8))
  abline(v=pars[par_no],lty=3, col=add.alpha('red',.85))
  xticks = pretty(range(c(lo,hi)))
  axis(side = 1, at = xticks*pars[par_no], xticks)
  mtext("Relative deviation from MLE", 1, line=2.25, cex=.75)
  if(par_no==1){
    ylbl = expression("l("~mu~"|"~gamma~","~beta~")")
  } else if(par_no==2){
    ylbl = expression("l("~gamma~"|"~mu~","~beta~")")
  } else if(par_no==3){
    ylbl = expression("l("~beta~"|"~mu~","~gamma~")")
  }  
  mtext(ylbl, 2, line=1, cex=.75)
  
}
ABC <- function(pars, Nt){
  
  u = pars[1]
  y = pars[2]
  b = pars[3]
  
  n = length(Nt)
  A = numeric(n)
  B = numeric(n)
  C = numeric(n)
  
  for(i in 2:n){
    ti    = Nt[i]-Nt[i-1]
    A[i] = exp( -b*ti )*( 1+A[i-1] )
    B[i] = ti*A[i]+exp(-b*ti)*B[i-1]
    C[i] = ti*(B[i]+exp(-b*ti)*B[i-1])+exp(-b*ti)*C[i-1]
  
    # tmp   = Nt[i]-Nt[1:(i-1)] # or  j    = which(Nt < Nt[i]); tmp  = Nt[i]-Nt[j]
    # A[i]  = sum( exp(-b*tmp) )
    # B[i]  = sum( exp(-b*tmp)*tmp )
    # C[i]  = sum( exp(-b*tmp)*tmp^2 )
  }
  return(data.frame(A=A,B=B,C=C))
}
get_se <- function(theta){
  H   = exp_hessian(pp,Tmax,theta)
  out = sqrt(diag(-solve(H)))
  return(out)
}
plot_SE <- function(par_no,pars,pp,Tmax, lo=.5, hi=1.5){
  # par_no = 1
  # pars = res
  # pp
  # Tmax
  # lo
  # hi
  
  nx   = 1000
  xmin = lo*pars[par_no]
  xmax = hi*pars[par_no]
  x = seq(xmin,xmax,length=nx)
  z = numeric(nx)
  for(i in 1:nx){
    if(par_no==1){
      z[i] = get_se(c(x[i],pars[2],pars[3]))[1]
    } else if(par_no==2){
      z[i] = get_se(c(pars[1],x[i],pars[3]))[2]
    } else if(par_no==3){
      z[i] = get_se(c(pars[1],pars[2],x[i]))[3]
    }
  }
  plot(x,z, xaxt='n', type='l', col=add.alpha('black',.75), bty='n', ann=FALSE)
  # abline(h=get_se(res)[par_no], col=add.alpha('black',.8))
  abline(v=pars[par_no],lty=3, col=add.alpha('red',.85))
  xticks = pretty(range(c(lo,hi)))
  axis(side = 1, at = xticks*pars[par_no], xticks)
  mtext("Relative deviation from MLE", 1, line=2.25, cex=.75)
  if(par_no==1){
    ylbl = expression("SE("~mu~"|"~gamma~","~beta~")")
  } else if(par_no==2){
    ylbl = expression("SE("~gamma~"|"~mu~","~beta~")")
  } else if(par_no==3){
    ylbl = expression("SE("~beta~"|"~mu~","~gamma~")")
  }  
  mtext(ylbl, 2, line=2.25, cex=.75)
}


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

# Make the gradients and Hessian of the likelihood:

# u: mu parameter (notation: \mu)
# y: branching ratio (notation: \gamma)
# b: beta parameter (notation: \beta)

(hgrad = hawkes_gradient(Tmax, pp, res))
exp_gradient(pp,Tmax,res)

lo = .75
hi = 1.25
par(mfrow=c(3,2), mar=c(4,3,1,1), oma=c(0,0,0,0))
plot_gradient(1,res,pp,Tmax,lo,hi) # Test gradient of mu
plot_SE(1,res,pp, Tmax, lo, hi)

plot_gradient(2,res,pp,Tmax,lo,hi) # Test gradient of gamma
plot_SE(2,res,pp, Tmax, lo, hi)

plot_gradient(3,res,pp,Tmax,lo,hi) # Test gradient of beta
plot_SE(2,res,pp, Tmax, lo, hi)




# Solve Hessian numerically
hhnum = hawkes_hessian(pp, Tmax, numerical=TRUE)
# Solve Hessian analytically
hhana = hawkes_hessian(pp, Tmax, numerical=FALSE)

ses.num = sqrt(diag(hhnum$I))
ses.ana = sqrt(diag(hhana$I))

# Compare analytical and numerical Hessian estimates:
ses.num/ses.ana
ses = ses.num
round(data.frame(mle = res, ci.lo = res-1.96*ses, ci.hi = res+1.96*ses, row.names = c("mu","gamma","beta")),4)
ses = ses.ana
round(data.frame(mle = res, ci.lo = res-1.96*ses, ci.hi = res+1.96*ses, row.names = c("mu","gamma","beta")),4)




