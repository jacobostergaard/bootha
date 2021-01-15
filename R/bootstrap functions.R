bootstrap_hawkes <- function(pp, pars, Tmax, B=299, type='fib', parametric=TRUE, dt=1, ncores=1){
  
  # Split the number of bootstrap samples on selected number of cores:
  b0 = ncores # desired number of batches
  b1 = B %% b0
  b2 = floor(B/b0)
  b  = c(rep(b2,b0)) # this is the list of smaller bootstrap samples to send to workers
  
  if(b1>0){ # not all B samples have been distributed in each batch, iterate until all have been placed
    btmp = b1
    i=1
    while(btmp>0){
      b[i] = b[i]+1 # place 1 sample in a batch and move on to the next one, until all are placed
      i = i+1
      btmp = btmp-1
    }
  }
  b = b[b>0]
  ncores = length(b)
  
  
  # Sampled waiting times (either parametric or non-parametric)
  
  
  if(type=='fib' | !parametric){
    x = seq(0,Tmax,dt)
    L = cbind(x,integrated_intensity(x, pp, pars))
    s = interpol(L, x=pp)
    v = c(s[1],diff(s)) 
    v = v/mean(v)
  } else{
    L = NULL
    s = NULL
    v = NULL
  }

  f <- function(i){
    
    out = lapply(1:i,function(i) {
                                    if(parametric){
                                      v_in = rexp(1.5*length(pp))  
                                    } else{
                                      v_in = sample(v,1.5*length(v),replace = TRUE)  
                                    }
                                # print(length(v_in))
                                    tmp = .bootstrap_aux(v_in, pars, Tmax, pp, type, L , dt)
                                # print(length(tmp))
                                    tmp = tmp[tmp<Tmax]
                                # print(length(tmp))
                                    return(tmp)
                                  })
    return(out)
  }
  
  if(ncores>1){ # Parallel part:
    cl = parallel::makeCluster(ncores, setup_timeout=.5)
      parallel::clusterExport(cl, c("pp","pars","Tmax","type","dt","L"), envir=environment())
      tmp = parallel::parLapply(cl, b, f)
    parallel::stopCluster(cl)  
    
    tmp = unlist(tmp, recursive = FALSE)
  
  } else{
    tmp = f(B)
  }
  
  return(tmp)
}


.bootstrap_aux <- function(v, pars, Tmax, pp, type='fib', L=NULL, dt=1){
  
  if(type=='fib'){
    if(is.null(L)){
      warning("Provide and integrated intensity for the Fixed Intensity Bootstrap")
    }
    # if(length(pars)==3){ # exponential kernel: use fast fib evaluation
      boot = exp_fib(L, v, pars, pp) # this is slower than interpolation...
    # } else{
      # boot = fib(L, v)
    # }
  } else{
    if(length(pars)==3){ # exponential kernel: use fast rib evaluation
      boot = exp_rib(pars, v)
    } else{
      boot = rib(pars, v, dt)  
    }
  }
  
  return(boot)
}
