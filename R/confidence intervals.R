boot.ci <- function(boot, pp, Tmax, pars, type='naive', a   = 0.05){
  # bootstrap confidence intervals
  
  # 'naive'   => find mle's for each bootstrap sample and find a/2, 1-a/2 quantiles
  # 'tratios' => find t-ratios based on the bootstrap sample mle's and Hessian and calculate the relevant quantiles
  # 'bootstd' => calculate the Hessian of each bootstrap sample, take the mean and use as std.error for asymptotic like CI's
  
  Nb = length(boot) # number of bootstrap samples
  
  f = function(pp){
    mle(pp, Tmax, init = pars, dt = 0)
  }
  
  g = function(pp){
    t_ratios(pp, Tmax, pars)
  }
  
  t_ratios <- function(pp, Tmax, pars){
    b.pars  = mle(pp,Tmax, pars)
    Hb      = exp_hessian(pp, Tmax, b.pars)
    Vb      = -solve(Hb)
    out     = (b.pars-pars)/sqrt(diag(Vb))
    return(out)
  }
  
  boot.mle      = matrix(unlist(lapply(boot, f)), nc=3, byrow=TRUE)
  pp.mle        = mle(pp, Tmax, init = pars, dt = 0)
  V             = -solve(exp_hessian(pp, Tmax, pp.mle))
  
  if(type=='naive'){
    qts = t(apply(boot.mle,2, function(x) quantile(x, probs=c(a/2, 1-a/2))))
    out = data.frame(mle=pp.mle, ci.lo = qts[,1], ci.hi = qts[,2])
    
  } else if(type=='tratios'){
    boot.t_ratios = matrix(unlist(lapply(boot, function(x){
                                                      b.pars  = mle(x,Tmax, pars)
                                                      Hb      = exp_hessian(x, Tmax, b.pars)
                                                      Vb      = -solve(Hb)
                                                      out     = (b.pars-pars)/sqrt(diag(Vb))
                                                      return(out)
                                                    })), nc=3, byrow = TRUE)
    
    qts    = t(apply(boot.t_ratios, 2, function(x) quantile(x,probs = c(1-a/2,a/2)))) # t ratio quantiles
    out    = data.frame(mle=pp.mle, ci.lo = pp.mle-sqrt(diag(V))*qts[,1], ci.hi = pp.mle-sqrt(diag(V))*qts[,2])
  } else if(type=='bootstd'){
    
    tmp = lapply(boot, function(x){
                            mle.tmp = mle(x, Tmax, pars, dt=0)
                            I = -solve(exp_hessian(x, Tmax, mle.tmp))
                      return(I)
                          } )
    VB  = array(unlist(tmp), dim=c(3,3,length(tmp)))
    VB  = apply(VB,c(1,2), mean)
    ses = sqrt(diag(VB))  
    
    qts = qnorm(1-a/2)
    out = data.frame(mle=pp.mle, ci.lo = pp.mle-qts*ses, ci.hi = pp.mle+qts*ses)
  }else{
    warning("Please specify type as either 'naive', 'tratios' or 'bootstd'")
  }
  
  return(out)
}

