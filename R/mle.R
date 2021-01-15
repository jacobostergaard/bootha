mle <- function(pp, Tmax, init, dt = 1, ui=NULL, ci=NULL, grad = NULL, hessian=FALSE, method = "Nelder-Mead"){
  
  if(length(init)==3){ # exponential kernel, use dt=0
    dt=0
  }
  
  f = function(init){
    -log_lik(pp, Tmax, init, dt)
  }
  if(is.null(ui)){
    ui = diag(length(init))
  }
  if(is.null(ci)){
    ci = rep(0,length(init))
  }
  if(!is.null(grad)){
    # If gradient is supplied, it is possible to obtain the Hessian
    mle = constrOptim(init, f, grad = grad, ui = ui, ci=ci, method=method, hessian=hessian)
    if(hessian){
      out = list(mle=mle$par, H=mle$hessian, I = -solve(mle$hessian) )  
    } else{
      out = mle$par
    }
  } else{
    # If no gradient is supplied, it is not possible to obtain the Hessian
    mle = constrOptim(init, f, ui = ui, ci=ci, method=method, hessian=FALSE)
    out = mle$par
  }
  
  # names(mle) = c("mu","n","B")
  return(out)
}

