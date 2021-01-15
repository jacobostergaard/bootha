kl_distance_chisq <- function(P,...){
  # Returns the Kullback-Leilber distance between an emperical distribution P and the chi^2(1) distribution
  
  h = hist(P, plot=FALSE, ...)  # Create a density at given midpoints to compare with theoretical distribution
  q = dchisq(h$mids,1)          # chi^(1) evaluated at midpoints from the density above
  p = h$density                 # the density from the histogram to evaluate against chi^2(1)
  
  # Fix any p=0, if there are areas with no observations in the empirical distribution, to avoid log(0) problem
  q = q[which(p>0)]   
  p = p[which(p>0)]

  out = sum(p*log(p/q))

  return(out)
}


kl_distance <- function(P,Q,...){
  # Returns the Kullback-Leilber distance between an emperical distributions P and Q
  
  # if(max(P)>max(Q)){
  #   h = hist(P, plot=FALSE, ...)              # Create a density at given midpoints to compare
  #   g = hist(Q, plot=FALSE, breaks=h$breaks)   
  # } else{
  #   h = hist(Q, plot=FALSE, ...)              # Create a density at given midpoints to compare
  #   g = hist(P, plot=FALSE, breaks=h$breaks)  
  # }
  
  h = hist(P, plot=FALSE, breaks=seq(-1,1.1*max(P,Q),length=500))              # Create a density at given midpoints to compare
  g = hist(Q, plot=FALSE, breaks=h$breaks)   
  
  p = h$density                 # the density from the histogram to evaluate against chi^2(1)
  q = g$density
  
  # Fix any p=0 or q=0, if there are areas with no observations in the empirical distribution, to avoid log(0), log(Inf) problem
  idx = which(p>0 & q>0)
  q   = q[idx]   
  p   = p[idx]
  
  out = sum(p*log(p/q))
  
  return(out)
}





plot_pp <- function(s,Tmax=NULL,...){
  if(is.null(Tmax)){
    Tmax = max(s)
  }
  # Plot a point process from input event times s and interval length [0,Tmax] 
  plot(c(0,s),1:(1+length(s)), type='s', ylim=c(0,length(s)), xlim=c(0,Tmax),...)
  points(s,rep(0,length(s)), pch=124)
}