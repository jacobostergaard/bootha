l.hawkes <- function(x,y, mu, a, b){
  
  out = numeric(length(x))  
  
  y   = sort(y)   # Make sure events are in order...
  nk  = length(y) # Events in input process
  
  # Evaluation points prior to any event
  idx = which(x<y[1])
  if(length(idx)>0){
    out[idx]  = mu  
  }
  
  # Evaluation points between first and last event
  for(i in 2:nk){
    idx = which(x<y[i] & x>=y[i-1])
    if(length(idx)>0){
      tmp = outer(x[idx],y[1:(i-1)],'-')
      tmp = apply(a*exp(-b*tmp),1,sum)
      out[idx] = mu+tmp
    }
  }
  
  # Evaluation points after last event
  idx = which(x>=y[nk])
  if(length(idx)>0){
    tmp = outer(x[idx],y,'-')
    tmp = apply(a*exp(-b*tmp),1,sum)
    out[idx] = mu+tmp  
  }
  
  if(is.null(y)){
    out = rep(mu,length(x))
  }
  
  return(out) 
}


L.hawkes <- function(x,y, mu, a, b){
  
  if(min(x)>0){
    dt    = mean(diff(x))
    x.tmp = c(seq(0,min(x)-dt,dt),x)
  } else{
    x.tmp = x
  }
  tmp = l.hawkes(x.tmp,y, mu, a, b)
  dx  = c(x[1],diff(x.tmp))
  out = cumsum(l.hawkes(x.tmp,y, mu, a, b))*dx
  
  idx = which(x.tmp >= min(x))
  out = out[idx]
  
  return(out)
}

# L.par <- function(x,y,a,b,mu,Lk){
#   if(y==0){ # First event
#     out = mu*x    
#   } else{ # Previous events noted
#     out = Lk+mu*(x-y)+a*(1-exp(-b*(x-y)))  
#   }
#   return(out)
# }
# L.aux <- function(x,y,Lk){
#   out = L.par(x,y,a,b,mu,Lk)
#   return(out)
# }
# 
# l.par <- function(x,y,mu,a,b,lk){
#   if(y==0){ # First event
#     out = rep(mu,length(x))
#   } else{
#     out = lk+a*exp(-b*(x-y))
#   }
#   return(out)
# }
# l.aux <- function(x,y,lk){
#   out = l.par(x,y,mu,a,b,lk)
# }
# 
# make_l <- function(pp, from, to, l, dt=0.001){
#   
#   idx = which(diff(pp)<dt) # Events with more than one in the dt size interval...
#   while(length(idx)>0){     
#     pp[idx+1] = pp[idx+1]+dt # add +dt to each of these to push into next interval
#     idx = which(diff(pp)<dt) # check if there are still problem intervals and continue untill this is not so
#   }
#   
#   
#   x     = seq(0,to,dt)      # evalutation points
#   i     = 1
#   idx   = x<=pp[i]          
#   l.tmp = l(x[idx],0,0)
#   l.tot = rep(NA,length(x))
#   l.tot[idx] = l.tmp
#   if(length(pp)>1){
#     for(i in 2:length(pp)){
#       z = l.tmp[length(l.tmp)]
#       idx = x>pp[i-1] & x<=pp[i]
#       l.tmp = l(x[idx],pp[i-1],z)
#       l.tot[idx] = l.tmp
#     }  
#   }
#   
#   z = l.tmp[length(l.tmp)]
#   idx = x>pp[i]
#   l.tmp = l(x[idx],pp[i],z)
#   l.tot[idx] = l.tmp  
#   
#   idx = which(x>=from)
#   out = data.frame(t = x[idx], l=l.tot[idx], L = (cumsum(l.tot)*dt)[idx])
#   
#   return(out)
# }
# 
# L = function(x,Tmax,dt=0.001){
#   tmp = make_l(hp,0,Tmax,l.aux,dt)
#   out = interp.L(L = data.frame(x=tmp$t,y=tmp$L), x=x)
#   return(out) 
# }
# l = function(x,Tmax,dt=0.001){
#   tmp = make_l(hp,0,Tmax,l.aux,dt)
#   out = interp.L(L = data.frame(x=tmp$t,y=tmp$l), x=x)
#   return(out) 
# }


