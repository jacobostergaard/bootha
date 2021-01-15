poisson_process <- function(Tmax,lam=1){
  # Returns a series of event times from a homogeneous Poisson process with parameter lam in the interval [0,T]
    
  
  # For future c++ implementation, the following is an example using U(0,1) variables rather than Exp(lam)
  # n = rpois(1,Tmax*lam)
  # U = runif(n)    # draw n U(0,1) random variables (for easier future c++ implementation)
  # E =-log(1-U)    # transform into Exp(1) variables
  # W = E/lam       # scale Exp(1) to Exp(lam) variables
  
  n = ceiling(Tmax*lam*2)   # approximate counts in [0,T] is lam*T, hence twice the expected counts, should avoid the while loop.
  W = rexp(n,lam)           # waiting times
  
  while(sum(W)<Tmax){       # If we still haven't reached the end of the interval [0,T] yet, simulate more waiting times
    W = rexp(n,lam) 
  }
  
  out = cumsum(W)   # event times
  out = out[out<Tmax]      # truncate to interval [0,T]
  
  return(out)
}


poisson_loglik <- function(n,Tmax,lam){
  # Returns the log-likelihood value of a Poisson process with n events in the interval [0,T] and parameter lam
  out = n*log(lam+1e-10)-lam*Tmax # the 1e-8 ensures a well defined log in the case of lam=0
  return(out)
}



logQ <- function(l1,l0,Tmax,n){
  # Returns the logQ, where Q is the ratio of two likelihoods with l1 and l0 parameters: log(L0/L1).
  l0 = poisson_loglik(n,Tmax,l0)
  l1 = poisson_loglik(n,Tmax,l1) 
  out  = l0-l1
  return(out)
}
