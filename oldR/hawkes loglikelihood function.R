llk <- function(theta, h, package="emhawkes", startval="mle"){
  # Return the log-likelihood value from input parameters and observed event times  
  if(h[1]>0){
    wobs = c(h[1],diff(h))  
  } else{
    wobs = diff(h)
  }
  
  if(startval=="mle"){
    llk1 = emhawkes::logLik(new("hspec", mu=theta[1], alpha=theta[2], beta=theta[3]), inter_arrival = wobs, lambda0 = theta[1])    
  } else{
    llk1 = misc::quiet(emhawkes::logLik(new("hspec", mu=theta[1], alpha=theta[2], beta=theta[3]), inter_arrival = wobs)  )
  }
  
  if(package == "emhawkes"){
    out = llk1  
  } else if(package == "hawkes"){
    llk2 = -hawkes::likelihoodHawkes(theta[1],theta[2],theta[3],h)
    out = llk2
  } else if(package == "both"){
    llk2 = -hawkes::likelihoodHawkes(theta[1],theta[2],theta[3],h)
    out = (llk1+llk2)/2
  }
  
  return(out)
}
