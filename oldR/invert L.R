# Invert integrated conditional intensity either as Fixed Intensity or Recursive Boostrap (FIB and RB)

invert.L <- function(x, L, Tmax, dt = 0.001){
  
  if(class(L)=="function"){
    # Recursive inversion
    # x = S.new
    # L = L.aux
    
    
    L.tmp = L(seq(0,Tmax,dt),0,0)
    tmp   = data.frame(x = seq(0,Tmax,dt), y = L.tmp-x[1])
    T.tmp = interp.L(tmp,y=0)
    T.out = T.tmp
    S.out = x[1]
    
    for(i in 2:length(x)){
      t1    = T.out[i-1]
      
      # t2    = seq(t1,Tmax,dt)
      if(t1<Tmax & !is.na(t1)){
        L.tmp = make_l(T.out,t1,Tmax, l = l.aux, dt = dt)  
        tmp   = data.frame(x = L.tmp$t, y = L.tmp$L-x[i])
        if(min(tmp$y)<0 & max(tmp$y>0)){
          T.tmp = interp.L(tmp,y=0)
          T.out = c(T.out,T.tmp)      
          S.out = c(S.out,x[i])
        }
      }
    }
    

    # # W     = x[1]
    # t2    = seq(0,Tmax,dt)
    # L.tmp = L(t2,0,0)
    # 
    # # S.tmp = W
    # S.tmp = x[1]
    # tmp   = data.frame(x = t2, y = L.tmp-S.tmp)
    # T.tmp = interp.L(tmp,y=0)
    # 
    # # T.out = T.tmp
    # # S.out = S.tmp
    # T.out = rep(T.tmp, length(x))
    # S.out = x
    # 
    # for(i in 2:length(x)){
    #   Lk    = S.tmp
    #   S.tmp = x[i]
    #   t1    = T.tmp
    #   t2    = seq(t1,Tmax,dt)
    #   L.tmp = L(t2,t1,Lk)
    #   tmp   = data.frame(x = t2, y = L.tmp-S.tmp)
    #   T.tmp = interp.L(tmp,y=0)
    #   
    #   T.out[i] = T.tmp
    # }
    # 
    
    # i = 2
    # while(T.tmp<Tmax){
    #   Lk    = S.tmp
    #   # W     = rexp(1)
    #   # S.tmp = S.tmp+W
    #   S.tmp = x[i]
    #   # S.out = c(S.out,S.tmp)
    #   t1    = T.tmp
    #   t2    = seq(t1,Tmax,dt)
    #   
    #   L.tmp = L(t2,t1,Lk)
    #   
    #   tmp   = data.frame(x = t2, y = L.tmp-S.tmp)
    #   if(max(tmp$y)>0){
    #     T.tmp = interp.L(tmp,y=0)
    #     T.out = c(T.out,T.tmp)
    #   } else{
    #     T.tmp = Tmax+1
    #   }
    #   i=i+1
    # }
    
  } else if(class(L)=="data.frame"){
    # Fixed intensity inversion
    S.out = x[x<=max(L[,2])] # Truncate input to observation interval
    T.out = interp.L(L ,y = S.out)
    
  }
  out   = data.frame(t=T.out, s=S.out) #[1:length(T.out)]) 

  return(out)
}


