bootstrap_sample <- function(pp,lb,Tmax){
  # Returns a bootstrapped point process from input process
  
  Nt  = length(pp)          # observed events 
  lh  = Nt/Tmax
  w   = c(pp[1],diff(pp))   # observed waiting times
  # v   = lb*w                # transformed waiting times using the bootstrap lambda
  if(is.null(lb)){
    v   = w
    lb  = 1
  } else{
    v   = w*lh  
  }
  
  vb  = sample(v, size=2*Nt, replace = TRUE)  # bootstrap sample in rescaled time (using twice the number of events in original process)
  
  while(sum(vb)<lb*Tmax){                     # If Tmax is not met, keep sampling until this is the case
    vb = c(vb,sample(v,1,replace=TRUE)) 
  }
  
  wb  = vb/lb           # bootstrap sample of waiting times
  
  ppb = cumsum(wb)      # bootstrapped point process
  ppb = ppb[ppb<= Tmax] # same, but truncated to [0,Tmax]
  
  return(ppb)  
}


bootstrap_pval <- function(pp,Tmax,l0=1,B=399){
  # Returns a p-value based on bootstrapped samples
  
  Nt = length(pp) # observed events 
  lh = Nt/Tmax    # MLE of lambda
  if(is.null(l0)){
    ppB  = lapply(rep(list(pp), B),bootstrap_sample,lb=NULL,Tmax=Tmax)  # B bootstrap samples of pp
  } else{
    ppB  = lapply(rep(list(pp), B),bootstrap_sample,lb=l0,Tmax=Tmax)  # B bootstrap samples of pp
  }
  
  NtB  = unlist(lapply(ppB, length))                                # Events in each sample
  
  # NtB  = NtB*l0/lh    # correction factor
  
  lhB  = NtB/Tmax                                                   # MLEs for each sample
  
  if(is.null(l0)){
    LRT  = -2*logQ(lh,lh,Tmax,Nt)   # LRT for observed process
    LRTB = -2*logQ(lhB,lh,Tmax,NtB) # LRTs for each sample
  } else{
    LRT  = -2*logQ(lh,l0,Tmax,Nt)   # LRT for observed process
    LRTB = -2*logQ(lhB,l0,Tmax,NtB) # LRTs for each sample
  }
  
  p_boot  = sum(LRTB>LRT)/B     # p-value bootstrap samples
  
  # hist(LRTB, breaks=100, prob=TRUE, border=NA, col='grey70')
  # curve(dchisq(x,1), add=TRUE, col='red', lwd=2)
  # abline(v=LRT, lwd=2)
  # 1-pchisq(LRT,1)
  
  return(p_boot)
}



bootstrap_pval_lamh <- function(pp,Tmax,l0,B=399,bootMLE=FALSE){
  # Returns a p-value based on bootstrapped samples
  
  Nt = length(pp) # observed events 
  lh = Nt/Tmax    # MLE of lambda
  
  if(bootMLE){
    lb = lh   # bootstrap using MLE of lambda
  } else{
    lb = l0   # bootstrap using true lambda
  }
  
  ppB  = lapply(rep(list(pp), B),bootstrap_sample,lb=lb,Tmax=Tmax)  # B bootstrap samples of pp
  NtB  = unlist(lapply(ppB, length))                                # Events in each sample
  
  # NtB  = NtB*l0/lh    # correction factor
  
  lhB  = NtB/Tmax                                                   # MLEs for each sample
  
  LRT  = -2*logQ(lh,l0,Tmax,Nt)   # LRT for observed process
  LRTB = -2*logQ(lhB,lh,Tmax,NtB) # LRTs for each sample
  
  p_boot  = sum(LRTB>LRT)/B     # p-value bootstrap samples
  
  # hist(LRTB, breaks=100, prob=TRUE, border=NA, col='grey70')
  # curve(dchisq(x,1), add=TRUE, col='red', lwd=2)
  # abline(v=LRT, lwd=2)
  # 1-pchisq(LRT,1)
  
  return(p_boot)
}
