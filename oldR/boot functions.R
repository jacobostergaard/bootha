# Bootstrap functions
pp.bootstrap1 <- function(L,s, kernel=NULL, quiet=FALSE, verbose=FALSE){
  L = as.matrix(L)
  
  if( !is.null(kernel) ){ # Recursive bootstrap
    out = NULL
    if(!quiet){
      print("Recursive bootstrap")  
    }
    for(i in 1:length(s)){
      f    = kernel(L[,1],out)-s[i]
      # tmp = suppressWarnings(interpol(cbind(L[,1],f),y=0))
      tmp  = interpol(cbind(L[,1],f),y=0)
      out  = c(out,tmp)
      if(verbose){ # Display progress
        cat("\r",round(100*i/length(s),2),"% done!     ")
      }
    }
    
  } else if(is.null(kernel)){ # Fixed intensity bootstrap
    if(!quiet){
      print("Fixed intensity bootstrap")  
    }
    out = interpol(L,y=s)
  }
  
  out = sort(out[!is.na(out)]) # Make sure event times are correctly ordered*
  return(out)  
  
  # *: this will almost surely be the case, however due to discretization of time + linear interpolation above, this may be necessary on rare occassions!
}


# pp.boot.old <- function(L,Smax,kernel,B, s.events = NULL, verbose=TRUE){
#   
#   
#   out = list()
#   for(i in 1:B){ # Parameteric bootstrap using exponential waiting times  
#     if(is.null(s.events)){
#       s     = poisson_process(Smax, 1) 
#     } else{ # Non-parametric bootstrap using input event times
#       
#       w   = c(s.events[1],diff(s.events)) # observed waiting times
#       W   = sample(w,size = 5*length(w), replace = TRUE) # sampled waiting times
#       s   = cumsum(W)  # Sampled oberved events
#   
#     }
#     
#     
#     b.tmp = pp.bootstrap1(L,s,kernel,TRUE,FALSE)  
#     
#     out[[i]] = b.tmp
#     if(verbose){ # Display progress
#       cat("\r",round(100*i/B,2),"% done!     ")
#     }
#   }
# 
#   return(out)
# }
# 
# 
# 
# pp.boot.par <- function(L,Smax,kernel,B, s.events = NULL, verbose=TRUE){
#   
#   library(parallel)
#   
#   f <- function(i){
#     if(is.null(s.events)){
#       s     = poisson_process(Smax, 1) 
#     } else{ # Non-parametric bootstrap using input event times
#       
#       w   = c(s.events[1],diff(s.events)) # observed waiting times
#       W   = sample(w,size = 5*length(w), replace = TRUE) # sampled waiting times
#       s   = cumsum(W)  # Sampled oberved events
#       
#     }
#     
#     b.tmp = pp.bootstrap1(L,s,kernel,TRUE,FALSE)
#     return(b.tmp)
#   }
#   
#   out = mclapply(1:B,f,mc.cores = 4)
# 
#   return(out)
# }


pp.boot <- function(L,Smax,kernel,B, s.events = NULL, parallel = FALSE, verbose=TRUE){
  
  
  out = list()
  
  f <- function(i=NULL){
    if(is.null(s.events)){
      s     = poisson_process(Smax, 1) 
    } else{ # Non-parametric bootstrap using input event times
      w   = c(s.events[1],diff(s.events)) # observed waiting times
      w   = w/mean(w) # scale waiting times so that the mean is 1
      
      W   = sample(w,size = 5*length(w), replace = TRUE) # sampled waiting times
      s   = cumsum(W)  # Sampled oberved events
    }
    
    out = pp.bootstrap1(L,s,kernel,TRUE,FALSE)
    out = out[out<=max(L[,1])]
    return(out)
  }
  
  if(parallel){
    out = parallel::mclapply(1:B,f,mc.cores = 4)
  } else{
    
    for(i in 1:B){ # Parameteric bootstrap using exponential waiting times  
      out[[i]] = f()
      if(verbose){ # Display progress
        cat("\r",round(100*i/B,2),"% done!     ")
      }
    }
    
  }
  
  return(out)
}


