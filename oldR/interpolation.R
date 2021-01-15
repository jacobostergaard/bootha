# Interpolation tools

interpolR <- function(f,x=NULL, y=NULL, flat=FALSE){
  # Interpolation of a two-dim function. 
  if(is.null(x)){
    fx = f[,2]
    fy = f[,1]
    
    x = y
  } else if(is.null(y)){
    fx = f[,1]
    fy = f[,2]
  }
  
  n = length(y)
  z = numeric(n)
  
  fx = fx[order(fx)]
  fy = fy[order(fx)]
  
  
  
  if(flat){ # Extrapolate extremes flat
    idx    = which(x <= min(fx))
    z[idx] = fy[1]
    idx    = which(x >= max(fx))
    z[idx] = fy[length(fy)]
  } else{ # Extrapolate extremes linearly
    N      = nrow(f)
    
    idx    = which(x <= min(fx))
    a      = diff(fy[1:2])/diff(fx[1:2])
    b      = fy[1]-fx[1]*a
    z[idx] = a*x[idx]+b
    
    idx    = which(x >= max(fx))
    a      = diff(fy[(N-1):N])/diff(fx[(N-1):N])
    b      = fy[N]-fx[N]*a
    z[idx] = a*x[idx]+b
  }
  
  # Interpolate linearly within extremes
  idx    = which(x > min(fx) & x < max(fx))
  for(i in idx){
    id.1 = max(which(fx<=x[i]))
    id.2 = min(which(fx>=x[i]))
    x.1  = fx[id.1]
    x.2  = fx[id.2]
    y.1  = fy[id.1]
    y.2  = fy[id.2]
    if(x.2!=x.1){
      w = (x.2-x[i])/(x.2-x.1)  
    } else{
      w = 1
    }
    z[i] = y.1*w+y.2*(1-w)
  }
  
  return(z)
}


interp.L <- function(L,x=NULL, y=NULL){
  # Linear interpolation of integrated conditoinal intensity, works both ways!
  
  if(is.null(x)){
    n = length(y)
    x = numeric(n)
    idx = which(y<=max(L[,2])) # Truncate y to interval of input function
    y   = y[idx]
    for(i in 1:n){
      id.1 = max(which(L[,2]<y[i]))
      id.2 = min(which(L[,2]>=y[i]))
      y.1 = L[id.1,2]
      y.2 = L[id.2,2]
      x.1 = L[id.1,1]
      x.2 = L[id.2,1]
      w = (y.2-y[i])/(y.2-y.1)
      x[i] = x.1*w+x.2*(1-w)
    }
    return(x)    
  } else if(is.null(y)){
    n = length(x)
    y = numeric(n)
    idx = which(x<=max(L[,1])) # Truncate x to interval of input function
    x   = x[idx]
    for(i in 1:n){
      id.1 = max(which(L[,1]<x[i]))
      id.2 = min(which(L[,1]>x[i]))
      x.1 = L[id.1,1]
      x.2 = L[id.2,1]
      y.1 = L[id.1,2]
      y.2 = L[id.2,2]
      w = (x.2-x[i])/(x.2-x.1)
      y[i] = y.1*w+y.2*(1-w)
    }
    return(y)    
  }
}