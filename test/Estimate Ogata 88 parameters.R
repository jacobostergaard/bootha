misc::clean_up()
library("PtProcess")

utsu = read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/UCPH/CIBS/Hawkes project/data/utsu.txt")

tmp = as.POSIXct(paste0(utsu$YEAR,"/",utsu$MO,"/",utsu$DY," ", utsu$HR,":",utsu$MN), tz="UTC")
time = as.numeric(difftime(tmp,tmp[1],units = "days"))
magnitude = utsu$MAG-6
 

# Shock types:
# 0 - main 
# 1 - fore
# 2 - after
ogata = data.frame(time,magnitude, shock = utsu$C) 
# ogata = ogata[ogata$shock==0,]
# ogata = ogata[ogata$shock!=1,]

# ogata$magnitude = ogata$magnitude*(ogata$shock==0)

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

ogata_gif <- function (data, evalpts, params, TT = NA, tplus = FALSE) {
  mu <- params[1]
  A <- params[2]
  alpha <- params[3]
  CC <- params[4]
  P <- params[5]
  if (is.null(data)) 
    data <- cbind(time = Inf, magnitude = 0)
  magnitudes <- data[, "magnitude"]
  times <- data[, "time"]
  if (any(is.na(TT))) {
    if (is.vector(evalpts)) 
      eval.times <- evalpts
    else eval.times <- evalpts[, "time"]
    if (length(eval.times) > 1) {
      St <- apply(t(eval.times), 2, function(z, magnitudes, 
                                             times, mu, A, P, CC, alpha, tplus) {
        if (!tplus) 
          use <- times < z
        else use <- times <= z
        if (sum(use) > 0) 
          sum(exp(alpha * magnitudes[use]) * (z - times[use]+CC)^(-P))
        else 0
      }, magnitudes = magnitudes, times = times, mu = mu, 
      A = A, P = P, CC = CC, alpha = alpha, tplus = tplus)
    }
    else {
      if (!tplus) 
        use <- times < eval.times
      else use <- times <= eval.times
      if (sum(use) > 0) 
        St <- sum(exp(alpha * magnitudes[use]) * (eval.times - times[use]+CC)^(-P))
      else St <- 0
    }
    ci <- mu + A * St
  }
  else {
    S <- c(NA, NA)
    for (Ii in 1:2) {
      Z0 <- (times < TT[Ii])
      if (any(Z0)) {
        if (P != 1) 
          S[Ii] <- sum(exp(alpha * data[, "magnitude"][Z0]) * 
                         (1 - (TT[Ii] - times[Z0]+CC)^(-P + 
                                                         1)))/(P - 1)
        else S[Ii] <- sum(exp(alpha * data[, "magnitude"][Z0]) * 
                            log(TT[Ii] - times[Z0] + CC))
      }
      else S[Ii] <- 0
    }
    ci <- mu * (TT[2] - TT[1]) + A * (S[2] - S[1])
  }
  names(ci) <- NULL
  return(ci)
}


TT <- c(0, max(ceiling(ogata$time)))

mu.in = nrow(ogata)/diff(range(TT))
alpha.in = 0.02
beta.in = 1.6
c.in =  0.02 # as.numeric(mean(diff(ogata$time)))
p.in = 1.01
params <- c(mu.in, alpha.in, beta.in, c.in, p.in, 1/mean(ogata$magnitude), 0) 


x <- mpp(data = ogata, gif = etas_gif, mark = list(dmagn_mark, rmagn_mark), params = params, TT = TT, gmap = expression(params[1:5]), mmap = expression(params))
# x <- mpp(data = ogata, gif = ogata_gif, mark = list(dmagn_mark, rmagn_mark), params = params, TT = TT, gmap = expression(params[1:5]), mmap = expression(params))

expmap <- function(y, p) {
  y$params[1:5] <- exp(p)
  return(y)
}
initial <- log(params[1:5])

hawkesboot::logQ
z <- optim(initial, neglogLik, object = x, pmap = expmap, control = list(trace = 1, maxit = 100))
 
initial <- z$par
z <- nlm(neglogLik, initial, object = x, pmap = expmap, print.level = 2, iterlim = 500, typsize = initial)
x0 <- expmap(x, z$estimate)

summary(x0)
print(logLik(x0))



res = residuals(x0)
plot(res,1:length(res), ylab = "Event Number", xlab = "Transformed Time", pty = "s", type='l')
abline(0,1, lty=3, col='red')


parest = as.data.frame(t(x0$params))
names(parest) = c("mu","K","beta","c","p","d1","d2")

parest$K = parest$K*parest$c^parest$p # Convert to Ogata par
parest = rbind(parest,c(0.00536, 0.017284,1.61385,  0.01959, 1.0, NA, NA))
rownames(parest) = c("Estimated", "Ogata 88")

round(parest,5)

