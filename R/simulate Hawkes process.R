simulateHawkes <- function(Tmax, pars, eps = 1e-10, verbose=FALSE){
  
  t = 0
  y = numeric(0)
  H = numeric(0)
  i = j = 1
  x = u = m = numeric(0)
  
  while(t < Tmax){
    M = intensity(t+eps, H, pars)
    E = rexp(1,M)
    t = t+E
    U = runif(1,0,M)
    
    if(verbose){
      x[j] = t
      u[j] = U
      m[j] = M
      j = j+1
      
    }
    
    if(t < Tmax & U < intensity(t,H,pars)){
      y[i] = t
      H = c(H,t)
      i=i+1
    }
  }
  
  if(!verbose){
    return(y)  
  } else{
    return(list(y=y, x=x, u=u, m=m))
  }
}


