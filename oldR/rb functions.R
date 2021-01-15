# Functions for the Recursive intensity bootstrap (RB)

misc::clean_up()

rb <- function(L, Tmax, dt = 0.001, type='par'){
  # Recursive bootstrap, currently only parametric bootstrap available
  # Here L is the integrated conditional intensity passed as a function call (parameters should be intrinsic to L)

  # Algorithm:
  #   1. Draw first waiting time Exp(1)
  #   2. Invert to oberserved time
  #   3. Recursively find the next oberserved waiting time by
  #       i)  Draw an Exp(1) variable
  #       ii) Invert the time interval using the integrated conditional intensity
  #   4. If the last observed time is after the input window:
  #       i) Stop and exit, removing the last waiting time
  #      Else go to 3.s

  W     = rexp(1)
  t2    = seq(0,Tmax,dt)
  L.tmp = L(t2,0,0)

  S.tmp = W
  tmp   = data.frame(x = t2, y = L.tmp-S.tmp)
  T.tmp = interp.L(tmp,y=0)

  T.out = T.tmp
  S.out = S.tmp

  while(T.tmp<Tmax){
    Lk    = S.tmp
    W     = rexp(1)
    S.tmp = S.tmp+W
    S.out = c(S.out,S.tmp)
    t1    = T.tmp
    t2    = seq(t1,Tmax,dt)

    L.tmp = L(t2,t1,Lk)

    tmp   = data.frame(x = t2, y = L.tmp-S.tmp)
    if(max(tmp$y)>0){
      T.tmp = interp.L(tmp,y=0)
      T.out = c(T.out,T.tmp)
      out   = data.frame(t=T.out, s=S.out)
    } else{
      T.tmp = Tmax+1
    }
  }

  return(out)

}

