# Functions to produce samples of waiting times in time transformed domain


# Parametric bootstrap: sample exponential waiting times
parboot_wt <- function(n, m=1){
  # Sample n waiting times, with mean m
  out = rexp(n,1/m)

  return(out)
}

resboot_wt <- function(in_wt, n = NULL, weights=NULL){
  # Resample input waiting times, n times (optional inputs) with weights (optional input)
  
  if(is.null(n)){
    n = length(in_wt)
  }
  if(is.null(weights)){
    out = sample(in_wt, size = n, replace = TRUE)
  } else{
    weights = weights/sum(weights) # Normalize weights
    out = sample(in_wt, size = n, replace = TRUE, prob = weights)
  }
  
}


