library(PtProcess)
library(emhawkes)
library(ppboot)
library(misc)
misc::clean_up()


mu   = .75
a    = .5
b    = 1.1
Tmax = 50
dt   = 1e-1

# Test error difference of fb and rb (they should be equal...)
set.seed(1234)
hp   = unlist(hawkes::simulateHawkes(mu,a,b,Tmax))
N    = length(hp)

x = seq(0,Tmax,dt)

mu.in = length(hp)/Tmax   # Set baseline guess to observed events/interval
a.in = mu.in/10           # Guess intensity jump is 10% of baseline
b.in = 2*a.in             # Ensure a/b ratio is less than 1

hspec <- new("hspec", mu=mu.in, alpha=a.in, beta=b.in)


hawkes_gif <- function(data, evalpts, params, TT=NA, tplus=FALSE){
  mu <- params[1]
  a <- params[2]
  b <- params[3]
  
  if (is.null(data)) 
    data <- cbind(time = Inf)
  times <- data[, "time"]
  
  if (any(is.na(TT))) {
    if (is.vector(evalpts)) 
      eval.times <- evalpts
    else eval.times <- evalpts[, "time"]
    if (length(eval.times) > 1) {
      St <- apply(t(eval.times), 2, function(z, times, mu, a, b, tplus) {
        if (!tplus) 
          use <- times < z
        else use <- times <= z
        if (sum(use) > 0) 
          sum(exp(-b * (z - times[use]) )) 
        else 0
      }, times = times, mu = mu, 
      a = a, b = b, tplus = tplus)
    }
    else {
      if (!tplus) 
        use <- times < eval.times
      else use <- times <= eval.times
      if (sum(use) > 0) 
        St <- sum(exp(-b * (eval.times - times[use]) ))
      else St <- 0
    }
    ci <- mu + a * St
  }
  else {
    S <- c(NA, NA)
    for (Ii in 1:2) {
      Z0 <- (times < TT[Ii])
      if (any(Z0)) {
        S[Ii] <- S[Ii] <- sum(exp(-b*(TT[Ii] - times[Z0])))
      }
      else S[Ii] <- 0
    }
    
    ci <- mu * (TT[2] - TT[1]) - a/b *(S[2] - S[1])
  }
  names(ci) <- NULL
  
  return(ci)  
}

dmagn_mark <- function(x, data, params) {
  if (params[7] > 0) {
    lambda <- etas_gif(data, x[, "time"], params = params[1:5])
    y <- dgamma(x[, "magnitude"], shape = 1 + sqrt(lambda) * params[7], rate = params[6], log = TRUE)
  } else y <- dexp(x[, "magnitude"], rate = params[6], log = TRUE) 
  
  return(y)
  
}

rmagn_mark <- function(ti, data, params) {
  if (params[7]>0) {
    lambda <- etas_gif(data, ti, params = params[1:5])
    y <- rgamma(1, shape = 1 + sqrt(lambda) * params[7], rate = params[6]) 
  } else y <- rexp(1, rate = params[6])
  
  return(list(magnitude = y))
}


invisible(capture.output({mle = hfit(hspec, inter_arrival = diff(hp), lambda0 = mu.in)}))  

mle$estimate




TT <- c(0, max(ceiling(hp)))

mu.in = mu
alpha.in = a
beta.in = b
params <- c(mu.in, alpha.in, beta.in) 

indat = data.frame(time = hp)

x <- mpp(data = indat, gif = hawkes_gif, mark = NULL, params = params, TT = TT, gmap = expression(params), mmap =NULL)

expmap <- function(y, p) {
  y$params <- exp(p)
  return(y)
}
initial <- log(params)

z <- optim(initial, neglogLik, object = x, pmap = expmap, control = list(trace = 1, maxit = 100))

initial <- z$par
z <- nlm(neglogLik, initial, object = x, pmap = expmap, print.level = 2, iterlim = 500, typsize = initial)
x0 <- expmap(x, z$estimate)

summary(x0)
print(logLik(x0))